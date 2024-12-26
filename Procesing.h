#pragma once
#include <array>
#include <iostream>
#include <ostream>
#include <cassert>
#include <random>
#include <tuple>
#include <ranges>
#include <cstring>
#include <algorithm>
#include "Fixed.h"
#include "RandomGenerator.h"

class ProcesingParent {
public:
    ProcesingParent() = default;
    ~ProcesingParent() = default;
    virtual void start(char** _field) = 0;
};


template<typename P, typename V, typename FLOW_TYPE, int N_SIZE, int M_SIZE>
class Procesing : public ProcesingParent{
public:
    static constexpr size_t N = N_SIZE;
    static constexpr size_t M = M_SIZE;
    static constexpr size_t T = 1'000'000;


    static constexpr std::array<std::pair<int,int>,4> deltas {{
        {-1,0},{1,0},{0,-1},{0,1}
    }};

    P rho[256];        
    P p[N][M]{};       
    P old_p[N][M];     
    int dirs[N][M]{};  
    char field[N][M + 1];

    template<typename VF>
    struct VectorField {
        std::array<VF,4> v[N][M];



        VF &add(int x,int y,int dx,int dy,VF dv){
            return get(x,y,dx,dy) += dv;
        }
        VF &get(int x, int y, int dx, int dy) {
            if(dx == -1 && dy == 0) { // Changed ranges to if
                return v[x][y][0];
            }
            else if(dx == 1 && dy == 0) {
                return v[x][y][1];
            }
            else if(dx == 0 && dy == -1) {
                return v[x][y][2];
            }
            else if(dx == 0 && dy == 1) {
                return v[x][y][3];
            }
            else {
                assert(false);
            }
        }
    };

    VectorField<V> velocity{};
    VectorField<FLOW_TYPE> velocity_flow{};

    

    int last_use[N][M]{};
    int UT = 0; 
        


        std::tuple<P, bool, std::pair<int, int>> propagate_flow(int x, int y, P lim) {
            last_use[x][y] = UT - 1;
            P ret = 0;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                    auto cap = velocity.get(x, y, dx, dy);
                    auto flow = velocity_flow.get(x, y, dx, dy);
                    if (flow == cap) {
                        continue;
                    }
                    // assert(v >= velocity_flow.get(x, y, dx, dy));
                    auto vp = std::min(lim, P(cap - flow));
                    if (last_use[nx][ny] == UT - 1) {
                        velocity_flow.add(x, y, dx, dy, vp);
                        last_use[x][y] = UT;
                        // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                        return {vp, 1, {nx, ny}};
                    }
                    auto [t, prop, end] = propagate_flow(nx, ny, vp);
                    ret += t;
                    if (prop) {
                        velocity_flow.add(x, y, dx, dy, t);
                        last_use[x][y] = UT;
                        // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                        return {t, prop && end != std::pair(x, y), end};
                    }
                }
            }
            last_use[x][y] = UT;
            return {ret, 0, {0, 0}};
        }

        void propagate_stop(int x, int y, bool force = false) {
            if (!force) {
                bool stop = true;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                        stop = false;
                        break;
                    }
                }
                if (!stop) {
                    return;
                }
            }
            last_use[x][y] = UT;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                    continue;
                }
                propagate_stop(nx, ny);
            }
        }

        P move_prob(int x, int y) {
            P sum = 0;
            for (size_t i = 0; i < 4; ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    continue;
                }
                sum += v;
            }
            return sum;
        }   

        struct ParticleParams {  
            char type;
            P cur_p;
            std::array<V, 4> v;

            void swap_with(Procesing &proc, int x, int y) {
                std::swap(proc.field[x][y], type);
                std::swap(proc.p[x][y], cur_p);
                std::swap(proc.velocity.v[x][y], v);
            }
        };

        bool propagate_move(int x, int y, bool is_first) {
            last_use[x][y] = UT - is_first;
            bool ret = false;
            int nx = -1, ny = -1;
            do {
                std::array<P, 4> tres;
                P sum = 0;
                for (size_t i = 0; i < 4; ++i) {
                    auto [dx, dy] = deltas[i];
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                        tres[i] = sum;
                        continue;
                    }
                    auto v = velocity.get(x, y, dx, dy);
                    if (v < 0) {
                        tres[i] = sum;
                        continue;
                    }
                    sum += v;
                    tres[i] = sum;
                }

                if (sum == 0) {
                    break;
                }

                P p = RandomGenerator::random01<P>() * sum;
                size_t d = tres.at(0) > p ? 0 : tres.at(1) > p ? 1 : tres.at(2) > p ? 2 :tres.at(3) > p ? 3 : 4; // Changed ranges to if

                auto [dx, dy] = deltas[d];
                nx = x + dx;
                ny = y + dy;
                assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

                ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
            } while (!ret);
            last_use[x][y] = UT;
            for (size_t i = 0; i < 4; ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                    propagate_stop(nx, ny);
                }
            }
            if (ret) {
                if (!is_first) {
                    ParticleParams pp{};
                    pp.swap_with(*this, x, y);
                    pp.swap_with(*this, nx, ny);
                    pp.swap_with(*this, x, y);
                }
            }
            return ret;
        }
    
        void start(char** _field) override {

            for (size_t x = 0; x < N; ++x)
            {   
                for (size_t y = 0; y < M; ++y)
                {
                    field[x][y] = _field[x][y];
                }
                field[x][M] = '\0';
            }

            rho[' '] = 0.01;
            rho['.'] = 1000;
            P g = 0.1;

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        dirs[x][y] += (field[x + dx][y + dy] != '#');
                    }
                }
            }

            for (size_t i = 0; i < T; ++i) {
                P total_delta_p = 0;
                // Apply external forces
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        if (field[x + 1][y] != '#')
                            velocity.add(x, y, 1, 0, g);
                    }
                }

                // Apply forces from p
                memcpy(old_p, p, sizeof(p));
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            int nx = x + dx, ny = y + dy;
                            if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                                auto delta_p = old_p[x][y] - old_p[nx][ny];
                                auto force = delta_p;
                                auto &contr = velocity.get(nx, ny, -dx, -dy);
                                if (contr * rho[(int) field[nx][ny]] >= force) {
                                    contr -= force / rho[(int) field[nx][ny]];
                                    continue;
                                }
                                force -= contr * rho[(int) field[nx][ny]];
                                contr = 0;
                                velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                                p[x][y] -= P(force) / P(dirs[x][y]);
                                total_delta_p -= P(force) / P(dirs[x][y]);
                            }
                        }
                    }
                }

                // Make flow from velocities
                velocity_flow = {};
                bool prop = false;
                do {
                    UT += 2;
                    prop = 0;
                    for (size_t x = 0; x < N; ++x) {
                        for (size_t y = 0; y < M; ++y) {
                            if (field[x][y] != '#' && last_use[x][y] != UT) {
                                auto [t, local_prop, _] = propagate_flow(x, y, 1);
                                if (t > 0) {
                                    prop = 1;
                                }
                            }
                        }
                    }
                } while (prop);

                // Recalculate p with kinetic energy
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] == '#')
                            continue;
                        for (auto [dx, dy] : deltas) {
                            auto old_v = velocity.get(x, y, dx, dy);
                            auto new_v = velocity_flow.get(x, y, dx, dy);
                            if (old_v > 0) {
                                assert(new_v <= old_v);
                                velocity.get(x, y, dx, dy) = new_v;
                                auto force = (old_v - new_v) * rho[(int) field[x][y]];
                                if (field[x][y] == '.')
                                    force *= 0.8;
                                if (field[x + dx][y + dy] == '#') {
                                    p[x][y] += force / P(dirs[x][y]);
                                    total_delta_p += force / P(dirs[x][y]);
                                } else {
                                    p[x + dx][y + dy] += P(force) / P(dirs[x + dx][y + dy]);
                                    total_delta_p += P(force) / P(dirs[x + dx][y + dy]);
                                }
                            }
                        }
                    }
                }

                UT += 2;
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            if (RandomGenerator::random01<P>() < move_prob(x, y)) {
                                prop = true;
                                propagate_move(x, y, true);
                            } else {
                                propagate_stop(x, y, true);
                               // std::cout << RandomGenerator::random01<P>() << " " << move_prob(x, y) << "\n";
                            }
                        }
                    }
                }

                if (prop) {
                    std::cout << "Tick " << i << ":\n";
                    for (size_t x = 0; x < N; ++x) {
                        std::cout << field[x] << "\n";
                    }
                }
            }
        }
};

