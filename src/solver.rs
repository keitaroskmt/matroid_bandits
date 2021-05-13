use std::collections::HashSet;
use rand_distr::{Normal, Distribution};
use rand::Rng;

use super::graph::*;

#[derive(Debug)]
pub struct GenerateReward {
    pub std_dev: f64,
}

impl GenerateReward {
    pub fn new() -> GenerateReward {
        GenerateReward {
            std_dev: 1.0,
        }
    }
    
    pub fn sample(&self, graph: &Vec<Vec<f64>>, e: &mut Edge, num: i64) -> f64 {
        let mean = graph[e.u][e.v]; 
        let normal = Normal::new(mean, self.std_dev).unwrap();
        let mut ret = 0.;
        for _ in 0..num {
            ret += normal.sample(&mut rand::thread_rng());
        }
        ret
    }
}

#[derive(Debug)]
pub struct Solver {
    graph: Vec<Vec<f64>>,
    es: Vec<Edge>,
    rank: i64,
    generate_reward: GenerateReward,
}

impl Solver {
    pub fn new(graph: Vec<Vec<f64>>, es: Vec<Edge>) -> Solver {
        let rank = Solver::calculate_rank(&es);
        let generate_reward = GenerateReward::new();
        let mut ret = Solver {
            graph,
            es,
            rank,
            generate_reward,
        };
        
        ret.init();
        ret
    }

    pub fn init(&mut self) {
        for e in &mut self.es {
            e.cost = 0.;
            e.sample_num = 0;
        }
    }
    
    pub fn sample_total(&self) -> i64 {
        let mut total = 0;
        for &e in &self.es {
            total += e.sample_num;
        }
        total
    }
    
    pub fn cost_diff(&self) -> Vec<f64> {
        let mut ret = Vec::new();
        for &e in &self.es {
            ret.push(e.cost - self.graph[e.u][e.v]);
        }
        ret
    }
    
    pub fn calculate_rank(es: &[Edge]) -> i64 {
        // calculate the number of connected components
        let mut vs = HashSet::new();
        for &e in es  {
            vs.insert(e.v);
            vs.insert(e.u);
        }
        let mut uf = UnionFind::new(vs.len());
        for &e in es {
            uf.unite(e.v, e.u);
        }
        let mut st = HashSet::new();
        for &i in &vs {
            st.insert(uf.root(i));
        }

        (vs.len() - st.len()) as i64
    }
    
    pub fn uniform_sample(&mut self, S: &HashSet<usize>, error_bound: f64, delta: f64) {
        for &i in S {
            let e = &mut self.es[i];
            let number = (error_bound.powi(-2) * (2.0 * delta.powi(-1)).ln() / 2.0) as i64;
            let value_sum = self.generate_reward.sample(&self.graph, e, number);
            e.cost = (e.sample_num as f64 * e.cost + value_sum) / (e.sample_num + number) as f64;
            e.sample_num += number;
        }
    }
    
    pub fn naive_i(&mut self, S: &HashSet<usize>, error_bound: f64, delta: f64) -> HashSet<usize> {
        self.uniform_sample(S, error_bound / 2.0, delta / (S.len() as f64));
        
        let mut es = Vec::new();
        for &i in S {
            es.push((i, self.es[i]));
        }
        optimal_base(self.es.len(), &mut es)
    }
    
    pub fn f_good_edges(&mut self, base_set: &HashSet<usize>, base: &HashSet<usize>, geta: f64) -> HashSet<usize> {
        // ret = {e \in base_set : base^{\geq mu_e + geta} does not block e}
        let mut ret = HashSet::new();

        // decreasing order
        let mut base_set_vec = base_set.into_iter().copied().collect::<Vec<usize>>();
        base_set_vec.sort_by(|&i, &j| self.es[j].cost.partial_cmp(&self.es[i].cost).unwrap());

        let mut base_vec = base.into_iter().copied().collect::<Vec<usize>>();
        base_vec.sort_by(|&i, &j| self.es[j].cost.partial_cmp(&self.es[i].cost).unwrap());

        let mut uf = UnionFind::new(self.es.len());
        let mut now = 0;
        for &i in &base_set_vec {
            while now < base.len() && self.es[base_vec[now]].cost >= self.es[i].cost + geta {
                uf.unite(self.es[base_vec[now]].u, self.es[base_vec[now]].v);
                now += 1;
            }
            if !uf.same(self.es[i].u, self.es[i].v) {
                ret.insert(i);
            }
        }
        
        ret
    }
    
    pub fn pac_sampleprune(&mut self, S: &HashSet<usize>, error_bound: f64, delta: f64) -> HashSet<usize> {
        let p: f64 = 0.01;
        if S.len() as f64 <= 2.0 * p.powi(-2) * (4.0 * (8.0 * delta.powi(-1)).ln()).max(self.rank as f64) {
            return self.naive_i(S, error_bound, delta)
        }
        
        // sample a subset F
        let mut rng = rand::thread_rng();
        let mut F = HashSet::new();
        for &i in S {
            if rng.gen_range(0.0..1.0) <= p {
                F.insert(i);
            }
        }

        let alpha = error_bound / 3.0;
        let lambda = error_bound / 12.0;
        let I = self.pac_sampleprune(&F, alpha, delta / 8.0);
        self.uniform_sample(S, lambda, delta * p / (8.0 * self.rank as f64));
        
        // s_next = I \cup {e \in S-I : I^{\geq mu_e - alpha - 2\lambda} does not block e}
        let s_minus_i = S.difference(&I).copied().collect::<HashSet<usize>>();
        let mut s_next = self.f_good_edges(&s_minus_i, &I, - alpha - 2.0 * lambda);
        s_next = s_next.union(&I).into_iter().copied().collect::<HashSet<usize>>();

        self.pac_sampleprune(&s_next, alpha, delta / 4.0)
    }
    
    pub fn exact_expgap(&mut self, S: &HashSet<usize>, error_bound: f64, delta: f64) -> HashSet<usize> {
        let mut r_elim = 1;
        let mut r_sele = 1;
        let mut s_cur = S.clone();
        let mut res = HashSet::new();
        loop {
            let mut es_cur = Vec::new();
            for &i in &s_cur {
                es_cur.push(self.es[i]);
            }
            let n_opt = Solver::calculate_rank(&es_cur);
            let n_bad = s_cur.len() as i64 - n_opt;
            println!("n_opt {:?}", n_opt);
            println!("n_bad {:?}", n_bad);
            
            if n_opt <= n_bad {
                if n_opt == 0 {
                    break;
                }
                let r = r_elim;
                let eps_r = 2.0f64.powi(-r) / 4.0;
                let delta_r = delta / (100.0 * (r as f64).powi(3));
                r_elim += 1;
                
                let I = self.pac_sampleprune(&s_cur, eps_r, delta_r);
                self.uniform_sample(&I, eps_r / 2.0, delta_r / n_opt as f64);
                let s_minus_i = s_cur.difference(&I).copied().collect::<HashSet<usize>>();
                self.uniform_sample(&s_minus_i, eps_r, delta_r / (n_opt as f64));
                
                let mut s_new = self.f_good_edges(&s_minus_i, &I, 1.5 * eps_r);
                s_new = s_new.union(&I).copied().collect::<HashSet<usize>>();
                s_cur = s_cur.intersection(&s_new).copied().collect::<HashSet<usize>>();

            } else {
                if n_bad == 0 {
                    res = res.union(&s_cur).copied().collect::<HashSet<usize>>();
                    break;
                }
                let r = r_sele;
                let eps_r = 2.0f64.powi(-r) / 4.0;
                let delta_r = delta / (100.0 * (r as f64).powi(3));
                r_sele += 1;

                self.uniform_sample(&s_cur, eps_r, delta_r / s_cur.len() as f64);
                
                // {e \in s_cur : (s_cur - {e})^{\geq mu_e - 2 * eps_r} does not block e}
                // cannot use f_good_edges because the base depends on the element of the base_set
                let mut U = HashSet::new();
                let mut s_cur_vec = s_cur.clone().into_iter().collect::<Vec<usize>>();
                s_cur_vec.sort_by(|&i, &j| self.es[j].cost.partial_cmp(&self.es[i].cost).unwrap());
                
                let mut uf = UnionFind::new(self.es.len());
                let mut now = 0;
                for &i in &s_cur_vec {
                    while now < s_cur.len() && self.es[s_cur_vec[now]].cost > self.es[i].cost {
                        uf.unite(self.es[s_cur_vec[now]].u, self.es[s_cur_vec[now]].v);
                        now += 1;
                    }
                    let uf_backup = uf.clone();
                    let now_backup = now;

                    // not contain e
                    now += 1;
                    while now < s_cur.len() && self.es[s_cur_vec[now]].cost > self.es[i].cost - 2.0 * eps_r {
                        uf.unite(self.es[s_cur_vec[now]].u, self.es[s_cur_vec[now]].v);
                        now += 1;
                    }

                    if !uf.same(self.es[i].u, self.es[i].v) {
                        U.insert(i);
                    }
                    
                    uf = uf_backup;
                    now = now_backup;
                    
                    res = res.union(&U).copied().collect::<HashSet<usize>>();
                    s_cur = self.matroid_contraction(&s_cur, &U);
                }
            }
        }
        res
    }
    
    pub fn matroid_contraction(&self, s: &HashSet<usize>, t: &HashSet<usize>) -> HashSet<usize> {
        let mut uf = UnionFind::new(self.es.len());
        for &e in t {
            uf.unite(self.es[e].u, self.es[e].v);
        }

        let mut ret = HashSet::new();
        for &e in s {
            if !uf.same(self.es[e].u, self.es[e].v) {
                ret.insert(e);
            }
        }
        
        ret
    }
    
    pub fn clucb(&mut self, S: &HashSet<usize>, delta: f64) -> HashSet<usize> {
        for &i in S {
            let e = &mut self.es[i];
            let val = self.generate_reward.sample(&self.graph, e, 1);
            e.cost = val;
            e.sample_num += 1;
        }
        
        let mut t = S.len() - 1;
        loop {
            t += 1;

            let mut es = Vec::new();
            for &i in S {
                es.push((i, self.es[i]));
            }
            let mt = optimal_base(S.len(), &mut es);
            
            let mut rads = vec![0.; S.len()];
            for &i in S {
                rads[i] = self.generate_reward.std_dev * (2.0 * (4.0 * S.len() as f64 * t as f64 / delta).ln() / self.es[i].sample_num as f64).sqrt();
                
                if mt.contains(&i) {
                    self.es[i].cost -= rads[i];
                } else {
                    self.es[i].cost += rads[i];
                }
            }
            println!("loop: {:?}, rad: {:?}", t, rads[0]);
            
            let mut es = Vec::new();
            for &i in S {
                es.push((i, self.es[i]));
            }
            let mt_new = optimal_base(self.es.len(), &mut es);
            
            if set_eq(&mt, &mt_new) {
                return mt;
            }
            
            
            let mut next_arm = 0;
            let mut mx = 0.;
            for &i in mt.symmetric_difference(&mt_new) {
                if mx > rads[i] {
                    mx = rads[i];
                    next_arm = i;
                }
            }
            
            let e = &mut self.es[next_arm];
            let val = self.generate_reward.sample(&self.graph, e, 1);
            e.cost = (e.cost * e.sample_num as f64 + val) / (e.sample_num + 1) as f64;
            e.sample_num += 1; 
        }
    }
}

fn set_eq(sa: &HashSet<usize>, sb: &HashSet<usize>) -> bool {
    (sa.len() == sb.len()) &&
        sa.iter().all(|i| sb.contains(i))
}
