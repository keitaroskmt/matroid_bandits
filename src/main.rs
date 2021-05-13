mod graph;
mod solver;

use std::collections::HashSet;
use rand_distr::{Normal, Distribution};
use rand::Rng;

use graph::*;
use solver::*;

// must be sorted
fn vec_eq(va: &[usize], vb: &[usize]) -> bool {
    (va.len() == vb.len()) &&
        va.iter().zip(vb).all(|(&a, &b)| a == b)
}


// give a large weight to only one spanning tree
// CLUCB can return optimal solution in one loop.
fn compare_sample1() {
    let n = 100;
    for &p in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        println!("\n------------");
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        for i in 0..n-1 {
            graph[i][i+1] = 10000.;
            graph[i+1][i] = 10000.;
            es.push(Edge::new(i, i+1, 1000.));
        }
        
        let mut rng = rand::thread_rng();
        for i in 0..n {
            for j in i+2..n {
                if rng.gen_range(0.0..1.0) <= p {
                    graph[i][j] = 10.;
                    graph[j][i] = 10.;
                    es.push(Edge::new(i, j, 10.));
                }
            }
        }
        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        
        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let error_bound = 1.;
        let delta = 0.1;
        
        // Chen, Lin, King, Lyu and Chen(2014)
        println!("-- CLUCB");
        let mut clucb_res = solver.clucb(&S, delta).into_iter().collect::<Vec<usize>>();
        clucb_res.sort();

        let sample_total = solver.sample_total();
        println!("sample_total {:?} when p = {:?}", sample_total, p);
        if vec_eq(&optimal, &clucb_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
        
        solver.init();
        // Chen, Gupta, and Li (2016)
        println!("\n-- Exact-ExpGap");
        let mut exact_res = solver.exact_expgap(&S, error_bound, delta).into_iter().collect::<Vec<usize>>();
        exact_res.sort();
        let sample_total = solver.sample_total();
        println!("sample_total {:?} when p = {:?}", sample_total, p);
        if vec_eq(&optimal, &exact_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
    }
}

// give randomized weight to edge
// Exact-ExpGap returns optimal solution quickly.
fn compare_sample2() {
    let n = 100;
    for &p in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        println!("\n------------");
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        let mut rng = rand::thread_rng();
        for i in 0..n {
            for j in i+1..n {
                if rng.gen_range(0.0..1.0) <= p {
                    let val = rng.gen_range(0.0..100000.0);
                    graph[i][j] = val;
                    graph[j][i] = val;
                    es.push(Edge::new(i, j, val));
                }
            }
        }

        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        
        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let error_bound = 1.;
        let delta = 0.1;
        
        // Chen, Lin, King, Lyu and Chen(2014)
        println!("-- CLUCB");
        let mut clucb_res = solver.clucb(&S, delta).into_iter().collect::<Vec<usize>>();
        clucb_res.sort();

        let sample_total = solver.sample_total();
        println!("sample_total {:?} when dencity = {:?}", sample_total, p);
        if vec_eq(&optimal, &clucb_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
        
        solver.init();
        // Chen, Gupta, and Li (2016)
        println!("\n-- Exact-ExpGap");
        let mut exact_res = solver.exact_expgap(&S, error_bound, delta).into_iter().collect::<Vec<usize>>();
        exact_res.sort();
        let sample_total = solver.sample_total();
        println!("sample_total {:?} when dencity = {:?}", sample_total, p);
        if vec_eq(&optimal, &exact_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
    }
}

fn main() {
    compare_sample1();
    compare_sample2();
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
        let mut es = Vec::new();
        es.push(Edge::new(0, 1, 2.));
        es.push(Edge::new(1, 2, 4.));
        es.push(Edge::new(0, 2, 3.));
        es.push(Edge::new(3, 4, 5.));
        let ret = Solver::calculate_rank(&es);
        println!("{}", ret);
        assert_eq!(ret, 3);
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
        
        solver.uniform_sample(&(0..n).map(|i| i).collect::<HashSet<usize>>(), 0.1, 0.1);
        println!("{:?}", solver);
    }

    #[test]
    fn test_sample_i() {
        // setting
        let n = 100;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        for i in 0..n-1 {
            graph[i][i+1] = 10000.;
            graph[i+1][i] = 10000.;
            es.push(Edge::new(i, i+1, 10000.));
        }
        
        for i in 0..n {
            for j in i+2..n {
                graph[i][j] = 10.;
                graph[j][i] = 10.;
                es.push(Edge::new(i, j, 10.));
            }
        }
        
        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        // [0, ..., n-2]
        assert_eq!(true, vec_eq(&optimal, &(0..n-1).collect::<Vec<usize>>()));


        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let error_bound = 1.;
        let delta = 0.1;
        
        println!("Naive-i");
        let mut collect_num = 0;

        let num = 10;
        for i in 0..num {
            solver.init();
            let mut naive_res = solver.naive_i(&S, error_bound, delta).into_iter().collect::<Vec<usize>>();
            naive_res.sort();
            
            let sample_total = solver.sample_total();
            println!("sample_total {:?}", sample_total);
            let cost_diff = solver.cost_diff();
            println!("max cost_diff {:?}", cost_diff.iter().fold(0./0., |s, t| t.max(s)));

            if vec_eq(&optimal, &naive_res) {
                collect_num += 1;
            }
        }
        println!{"error rate {:?}", 1.0 - collect_num as f64 / num as f64};
    }

    #[test]
    fn test_pac_sampleprune() {
        // setting
        let n = 100;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        for i in 0..n-1 {
            graph[i][i+1] = 10000.;
            graph[i+1][i] = 10000.;
            es.push(Edge::new(i, i+1, 10000.));
        }
        
        for i in 0..n {
            for j in i+2..n {
                graph[i][j] = 10.;
                graph[j][i] = 10.;
                es.push(Edge::new(i, j, 10.));
            }
        }
        
        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        // [0, ..., n-2]
        assert_eq!(true, vec_eq(&optimal, &(0..n-1).collect::<Vec<usize>>()));
        

        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let error_bound = 1.;
        let delta = 0.1;

        println!("PAC-SamplePrune");

        let mut pac_res = solver.pac_sampleprune(&S, error_bound, delta).into_iter().collect::<Vec<usize>>();
        pac_res.sort();
        
        let sample_total = solver.sample_total();
        println!("sample_total {:?}", sample_total);
        let cost_diff = solver.cost_diff();
        println!("max cost_diff {:?}", cost_diff.iter().fold(0./0., |s, t| t.max(s)));

        if vec_eq(&optimal, &pac_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
    }
    
    #[test]
    fn test_exact_expgap() {
        // setting
        let n = 100;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        for i in 0..n-1 {
            graph[i][i+1] = 10000.;
            graph[i+1][i] = 10000.;
            es.push(Edge::new(i, i+1, 10000.));
        }
        
        for i in 0..n {
            for j in i+2..n {
                graph[i][j] = 10.;
                graph[j][i] = 10.;
                es.push(Edge::new(i, j, 10.));
            }
        }
        
        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        // [0, ..., n-2]
        assert_eq!(true, vec_eq(&optimal, &(0..n-1).collect::<Vec<usize>>()));
        

        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let error_bound = 1.;
        let delta = 0.1;

        println!("Exact-ExpGap");

        let mut exact_res = solver.exact_expgap(&S, error_bound, delta).into_iter().collect::<Vec<usize>>();
        exact_res.sort();
        
        let sample_total = solver.sample_total();
        println!("sample_total {:?}", sample_total);
        let cost_diff = solver.cost_diff();
        println!("max cost_diff {:?}", cost_diff.iter().fold(0./0., |s, t| t.max(s)));

        if vec_eq(&optimal, &exact_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
    }

    #[test]
    fn test_clucb() {
        // setting
        let n = 100;
        let mut graph = vec![vec![0.; n]; n];
        let mut es = Vec::new();
        let mut rng = rand::thread_rng();
        for i in 0..n {
            for j in i+1..n {
                if rng.gen_range(0.0..1.0) <= 0.8 {
                    let val = rng.gen_range(0.0..100000.0);
                    graph[i][j] = val;
                    graph[j][i] = val;
                    es.push(Edge::new(i, j, val));
                }
            }
        }
        
        // optimal
        let mut es_idx = Vec::new();
        for (i, e) in es.iter().enumerate() {
            es_idx.push((i, e.clone()));
        }
        let mut optimal = optimal_base(n, &mut es_idx).into_iter().collect::<Vec<usize>>();
        optimal.sort();
        

        let size = es.len();
        let mut solver = Solver::new(graph, es);
        let S = (0..size).collect::<HashSet<usize>>();
        let delta = 0.1;

        println!("CLUCB");

        let mut clucb_res = solver.clucb(&S,delta).into_iter().collect::<Vec<usize>>();
        clucb_res.sort();
        
        let sample_total = solver.sample_total();
        println!("sample_total {:?}", sample_total);
        let cost_diff = solver.cost_diff();
        println!("max cost_diff {:?}", cost_diff.iter().fold(0./0., |s, t| t.max(s)));

        if vec_eq(&optimal, &clucb_res) {
            println!("Correct")
        } else {
            println!("Wrong")
        }
    }
}