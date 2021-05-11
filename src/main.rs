mod graph;
mod solver;

use std::collections::HashSet;
use rand_distr::{Normal, Distribution};
use rand::Rng;

use graph::*;
use solver::*;


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
}