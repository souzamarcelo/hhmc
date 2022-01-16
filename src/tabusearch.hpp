#pragma once

struct BTR {
    double p;
    unsigned i;
    
    BTR(double _p) : p(_p), i(0) {}
    
    // some improvement w/ reservoir sampling
    template<typename Excluder>
    unsigned sie(const Solution &S, Excluder exclude, int bkv) {
        unsigned imp = S.x.size(), emp = S.x.size();
        unsigned nimp = 0, nemp = 0;
        int nt = (S.x.size() + 9) / 10;
        unsigned tt = S.x.size();
        while (nt-- > 0 || (tt-- > 0 && nimp == 0 && nemp == 0)) {
            // aspiration?
            if(S.flipvalue(i) + S.value < bkv)
                return i;
            if(!exclude(i) && S.flipvalue(i) < 0) {
                nimp++;
                if(nimp * random01(rng) <= 1)
                    imp = i;
            }
            if(!exclude(i) && S.flipvalue(i) <= 0) {
                nemp++;
                if(nemp * random01(rng) <= 1)
                    emp = i;
            }
            i = (i + 1) % S.x.size();
        }
        if(nimp > 0)
            return imp;
        else return emp;
    }
    
    template<typename Excluder>
    unsigned step(const Solution &S, Excluder exclude, int bkv) {
        if(random01(rng) < p)
            return randomInt(rng) % S.x.size();
        else {
            unsigned r = sie(S, exclude, bkv);
            if(isValid(r, S))
                return r;
            else
                return randomInt(rng) % S.x.size();
        }
    }
};

template<typename Improvement, typename Duration>
unsigned tabusearch(Solution &S, Improvement improve, Duration dgen, chrono::system_clock::time_point start, int target,
                    const unsigned maxstagnate = numeric_limits<unsigned>::max(),
                    const unsigned maxsteps = numeric_limits<unsigned>::max()) {

    unsigned i, steps = 0, notimproved = 0;
    Solution B(S.I);
    B = S;

    report.newBestKnownValue(chrono::system_clock::now(), B.value);

    vector<unsigned> tabu(S.x.size() + 1, 0);
    
    do {
        if(S.value <= target)
            break;
        if(chrono::system_clock::now() - start > chrono::duration<int>(timeLimit))
            break;
        if(notimproved > maxstagnate)
            break;
        if(steps > maxsteps)
            break;
        
        steps++;
        i = improve.step(S, [&steps, &tabu](unsigned j) { return steps < tabu[j]; }, B.value);
        if(!isValid(i, S))
            continue;
        S.flip(i);
        
        unsigned d = dgen();
        if(d == S.x.size())
            tabu[i] = steps + S.I.deg[i] + 1;
        else if(d < S.x.size())
            tabu[i] = steps + d + 1;
        
        if(S.value < B.value) {            
            B = S;
            B.time = chrono::system_clock::now();
            notimproved = 0;

            report.newBestKnownValue(chrono::system_clock::now(), B.value);
        } else
            notimproved++;
    } while (true);
    steps--;
    S = B;
    return steps;
}
