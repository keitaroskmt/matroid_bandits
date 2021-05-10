mod kruskal;

use rand_distr::{Normal, Distribution};
use rand::Rng;
use kruskal::*;

#[derive(Debug)]
pub struct GenerateReward {
    std_dev: f64,
}

impl GenerateReward {
    pub fn new() -> GenerateReward {
        GenerateReward {
            std_dev: 1.0,
        }
    }
    
    pub fn sample(&self, graph: &Vec<Vec<f64>>, e: &mut Edge) -> f64 {
        let mean = graph[e.u][e.v]; 
        let normal = Normal::new(mean, self.std_dev).unwrap();
        let ret = normal.sample(&mut rand::thread_rng());
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
        let rank = Solver::calculate_rank(&graph);
        let generate_reward = GenerateReward::new();
        Solver {
            graph,
            es,
            rank,
            generate_reward,
        }
    }
    
    fn calculate_rank(graph: &[Vec<f64>]) -> i64 {
        // calculate the number of connected components
        let mut seen = vec![false; graph.len()];
        let mut ret = 0;
        for i in 0..(graph.len()) {
            if !seen[i] {
                Solver::dfs(graph, &mut seen, i);
                ret += 1
            }
        }
        ret = graph.len() as i64 - ret;
        ret
    }
    
    fn dfs(graph: &[Vec<f64>], seen: &mut Vec<bool>, v: usize) {
        seen[v] = true;
        // next vertex
        let nvs = graph[v].iter()
                    .enumerate()
                    .filter(|&(_, &x)| x > 0.0)
                    .map(|(i, _)| i);
    
        for nv in nvs {
            if !seen[nv] {
                Solver::dfs(graph, seen, nv);
            }
        }
    }
    
    fn uniform_sample(&mut self, idx: &[usize], error_bound: f64, delta: f64) {
        for &i in idx {
            let e = &mut self.es[i];
            let number = (error_bound.powi(-2) * (2.0 * delta.powi(-1)).ln() / 2.0) as i64;
            for _ in 0..number {
                let value = self.generate_reward.sample(&self.graph, e);
                e.cost = (e.sample_num as f64 * e.cost + value) / (e.sample_num + 1) as f64;
                e.sample_num += 1;
            }
        }
    }
    
    fn naive_i(&mut self, idx: &[usize], error_bound: f64, delta: f64) -> Vec<usize> {
        self.uniform_sample(idx, error_bound / 2.0, delta / (idx.len() as f64));
        
        let mut es = Vec::new();
        for &i in idx {
            es.push((i, self.es[i]));
        }
        optimal_base(es.len(), &mut es)
    }
    
    fn pac_sampleprune(&mut self, idx: &[usize], error_bound: f64, delta: f64) -> Vec<usize> {
        let p: f64 = 0.01;
        if idx.len() as f64 <= 2.0 * p.powi(-2) * (4.0 * (8.0 * delta.powi(-1)).ln()).max(self.rank as f64) {
            return self.naive_i(idx, error_bound, delta)
        }
        
        // sample a subset F
        let mut rng = rand::thread_rng();
        let mut F = Vec::new();
        for &i in idx {
            if rng.gen_range(0.0..1.0) <= p {
                F.push(i)
            }
        }

        let alpha = error_bound / 3.0;
        let lambda = error_bound / 12.0;
        let I = self.pac_sampleprune(&F, alpha, delta / 8.0);
        self.uniform_sample(idx, lambda, delta * p / (8.0 * self.rank as f64));
        todo!()
    }
}


fn solve() {
}

fn main() {
    solve()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn normal() {
        let normal = Normal::new(2.0, 3.0).unwrap();
        let v = normal.sample(&mut rand::thread_rng());
        println!{"{} is from N(2, 9) distribution", v};
    }

    #[test]
    fn es_sort() {
        let mut es = Vec::new();
        es.push(Edge::new(0, 1, 2.));
        es.push(Edge::new(1, 2, 4.));
        es.push(Edge::new(0, 2, 3.));
        println!("{:?}", es);
        es.sort_by(|s, t| t.cost.partial_cmp(&s.cost).unwrap()); 
        println!("{:?}", es);
    }

    #[test]
    fn test_rank() {
        let n = 5;
        let mut graph = vec![vec![0.; n]; n];
        graph[0][1] = 1.;
        graph[1][2] = 2.;
        graph[2][0] = 3.;
        graph[3][4] = 3.;
        let ret = Solver::calculate_rank(&graph);
        println!("{}", ret);
        assert_eq!(ret, 2);
    }

    #[test]
    fn test_kruskal() {
        let n = 3;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        
        graph[0][1] = 5.;
        graph[1][0] = 5.;
        es.push(Edge::new(0, 1, 5.));

        graph[1][2] = 10.;
        graph[2][1] = 10.;
        es.push(Edge::new(1, 2, 10.));

        graph[2][0] = 1.;
        graph[0][2] = 1.;
        es.push(Edge::new(0, 2, 1.));

        let res = kruskal(n, &mut es);
        assert_eq!(res, 6.);
    }

    #[test]
    fn test_struct() {
        let n = 3;
        let mut es = Vec::new();
        es.push(Box::new(Edge::new(0, 1, 5.)));
        es.push(Box::new(Edge::new(1, 2, 10.)));
        es.push(Box::new(Edge::new(0, 2, 1.)));
        for e in &mut es {
            e.sample_num += 1;
            println!("{:?}", e);
        }

        // sample a subset F
        let mut F = Vec::new();
        for e in &es {
            F.push(e.clone())
        }
        for e in &mut F {
            e.sample_num += 1;
            println!("{:?}", e);
        }

        for e in &mut es {
            println!("{:?}", e);
        }
    }

    #[test]
    fn test_sample() {
        let n = 3;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        
        graph[0][1] = 5.;
        graph[1][0] = 5.;
        es.push(Edge::new(0, 1, 5.));

        graph[1][2] = 10.;
        graph[2][1] = 10.;
        es.push(Edge::new(1, 2, 10.));

        graph[2][0] = 1.;
        graph[0][2] = 1.;
        es.push(Edge::new(0, 2, 1.));
        

        let mut solver = Solver::new(graph, es);
        println!("{:?}", solver);
        
        solver.uniform_sample(&(0..n).map(|i| i).collect::<Vec<usize>>(), 0.1, 0.1);
        println!("{:?}", solver);
    }
}