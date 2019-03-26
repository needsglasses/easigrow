//! Numerical optimisation using the Nelder-Mead algorithim
//!
//! Translation of the pure Python/Numpy implementation of the Nelder-Mead algorithm.
//! Reference: <https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method>

// We use this crate to allow us to add and multiply vectors to keep it in line with python.
// similar to nalgebra but nalgrebra requires way too many subcrates
use vector::MVector;
use log::{info, warn};

/// A point on the simplex.
#[derive(Debug, Clone)]
struct Result {
    x: MVector<f64>,
    score: f64,
}

// finds the centroids of a list of points
fn centroid(simplex: &[Result]) -> MVector<f64> {
    let n = simplex[0].x.len();
    let mut zero = Vec::new();

    for _ in 0..n {
        zero.push(0.0f64);
    }

    //    let x0 = MVector::from_slice(n, &zero);
    let x0 = MVector::from_row_slice(n, &zero);
    simplex.iter().fold(x0, |sum, r| sum + r.x.clone()) / simplex.len() as f64
}

// Keep track of each type of modification.
#[derive(Debug)]
enum Operation {
    Reflection,
    Expansion,
    Contraction,
    Reduction,
}

/// Nelder-Mead specific optimisation parameters.
#[derive(Debug, Clone)]
pub struct Nelder {
    /// step increment in each direction from the starting point
    pub step: f64,
    /// reflection factor (relects the worst point through the centroid)
    pub alpha: f64,
    /// expansion factor (extends the reflected point)
    pub gamma: f64,
    /// contraction factor (away from the worst point)
    pub rho: f64,
    /// shrinking factor (around the best point)
    pub sigma: f64,
}

impl Nelder {
    pub fn new(x: &[f64]) -> Nelder {
        Nelder {
            step: x[0],
            alpha: x[1],
            gamma: x[2],
            rho: x[3],
            sigma: x[4],
        }
    }
    pub fn default() -> Nelder {
        Nelder {
            step: 0.1,
            alpha: 1.0,
            gamma: 2.0,
            rho: 0.2,
            sigma: 0.5,
        }
    }
}

//// Nelder-Mead non-linear optimisation
///
/// This routine works best if each of the parameters being optimised are
/// roughly the same size.  If this is not the case then they should
/// be normalised to ensure they are.
pub fn nelder_mead<F>(
    f: &F,
    x: &mut [f64],
    params: &Nelder,
    converge_tol: f64,
    max_iter: usize,
) -> f64
where
    F: Fn(&[f64]) -> f64,
{
    check_nelder_limits(params);

    let x_start = MVector::from_row_slice(x.len(), x);
    let mut results = Vec::new();
    let mut ops = Vec::new();

    // results.push(Result {x: x_start.clone(), score: f(x_start.as_ref())});
    results.push(Result {
        x: x_start.clone(),
        score: f(x_start.as_slice()),
    });

    // iniitalise the simplex by taking a step in each direction.
    for i in 0..x_start.len() {
        let mut x_init = x_start.clone();
        // just have to peer inside vec here cause I can't make it work otherwise
        x_init.m[i] += params.step;
        // results.push(Result {x: x_init.clone(), score: f(x_init.as_ref())});
        results.push(Result {
            x: x_init.clone(),
            score: f(x_init.as_slice()),
        });
    }
    info!("Nelder: starting result {:?}", results);

    let mut prev_best = results[0].score;
    let mut iter = 0;
    let n = results.len();

    loop {
        // Check if exceeding the iteration limit.
        if iter >= max_iter {
            warn!("***Warning: The optimisation has failed to converge within the specified maximum iteration limit {}. 
The answer may not be optimum. Try increasing the limit.", max_iter);
            break;
        }
        iter += 1;

        results.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
        info!("Nelder: {} {}. best: {}", iter, n, results[0].score);

        // check for convergence
        let change_tol = (results[0].score - results[n - 1].score).abs();
        if change_tol < converge_tol {
            info!("Nelder: Success. Converged tol {}", change_tol);
            break;
        } else {
            info!("Nelder: convergence tolerance not reached {}", change_tol);
            info!(
                "Nelder: best {}, worst {}, prev_best {}",
                results[0].score,
                results[results.len() - 1].score,
                prev_best
            );
        }

        prev_best = results[0].score;

        // calculate centro
        let x0 = centroid(&results[0..n - 1]);

        // reflection
        // if the reflected point is better than the second worst,
        // but not better than the best, then replace the worst point
        // with the reflected point
        let xr = x0.clone() + params.alpha * (x0.clone() - results[n - 1].x.clone());
        // let rscore = f(xr.as_ref());
        let rscore = f(xr.as_slice());
        info!("Nelder: rscore: {}", rscore);

        if rscore >= results[0].score && results[n - 2].score > rscore {
            ops.push(Operation::Reflection);
            info!("Nelder: Including reflected point");

            results.pop().unwrap();
            results.push(Result {
                x: xr,
                score: rscore,
            });
            continue;
        }

        // expansion
        // If the reflected point is the best point so far
        // then compute the expanded point
        if rscore < results[0].score {
            let xe = x0.clone() + params.gamma * (xr.clone() - x0.clone());
            // let escore = f(xe.as_ref());
            let escore = f(xe.as_slice());
            info!("Nelder: expansion score: {}", escore);

            results.pop().unwrap();

            if escore < rscore {
                ops.push(Operation::Expansion);
                results.push(Result {
                    x: xe,
                    score: escore,
                });
                continue;
            } else {
                results.push(Result {
                    x: xr,
                    score: rscore,
                });
                continue;
            }
        }

        // contraction
        // If the contracted point is better than the worst point,
        // replace the worst point with the contracted point
        let xc = x0.clone() + params.rho * (results[n - 1].x.clone() - x0);
        // let cscore = f(xc.as_ref());
        let cscore = f(xc.as_slice());

        if cscore < results[n - 1].score {
            info!("contracting: {}", cscore);

            ops.push(Operation::Contraction);
            results.pop().unwrap();
            results.push(Result {
                x: xc,
                score: cscore,
            });
            continue;
        }

        // reduction
        // For all but the best point, replace the point with
        // xi = x1 + sigma(xi - x)
        // This is a shrinking of the simplex around the best point
        info!("Nelder: reducing");

        ops.push(Operation::Reduction);
        for r in 1..results.len() {
            results[r].x =
                results[0].x.clone() - params.sigma * (results[r].x.clone() - results[0].x.clone());
            // results[r].score = f(results[r].x.as_ref());
            results[r].score = f(results[r].x.as_slice());
        }
    }

    info!("Nelder: Iterations: {}", iter);
    let count = count_nelder_ops(&ops);
    info!(
        "Nelder: Operations: reflection {}, expansion {}, contraction {}, reduction {}",
        count.reflection, count.expansion, count.contraction, count.reduction
    );
    for (i, &a) in results[0].x.clone().as_slice().to_vec().iter().enumerate() {
        x[i] = a;
    }
    info!("nelder x: {:?}", x);
    results[0].score
}

struct NelderSum {
    reflection: usize,
    expansion: usize,
    contraction: usize,
    reduction: usize,
}

// Count the number of each type of operation in the nelder search.
fn count_nelder_ops(ops: &[Operation]) -> NelderSum {
    ops.iter().fold(
        NelderSum {
            reflection: 0,
            expansion: 0,
            contraction: 0,
            reduction: 0,
        },
        |mut count, op| {
            match *op {
                Operation::Reflection => count.reflection += 1,
                Operation::Expansion => count.expansion += 1,
                Operation::Contraction => count.contraction += 1,
                Operation::Reduction => count.reduction += 1,
            };

            count
        },
    )
}

/// Check that the Nelder parameters are acceptable.
fn check_nelder_limits(params: &Nelder) {
    // recommended parameter limits (wikipedia)
    if params.gamma < 0.0 {
        warn!("***Warning: Wikipedia recommends using a Nelder-Mead value for gamma > 0.0, using {} .", params.gamma);
    }
    if { params.rho < 0.0 } | { params.rho > 0.5 } {
        warn!(
            "***Warning: Wikipedia recommends using a Nelder-Mead value 0.0 < rho < 0.5, using {}.",
            params.rho
        );
    }
    if params.sigma < 0.0 {
        warn!(
            "***Warning: Wikipedia recommends using a Nelder-Mead value for sigma > 0.0, using {}.",
            params.sigma
        );
    }
}

#[cfg(test)]
mod tests {
    use vector::MVector;
    use log::info;

    use super::{centroid, nelder_mead};
    use super::{Nelder, Result};

    #[test]
    fn test_centroid() {
        let x = vec![
            Result {
                x: MVector::from_row_slice(3, &[0.0, 0.0, 0.0]),
                score: 0.0,
            },
            Result {
                x: MVector::from_row_slice(3, &[1.0, 2.5, 1.0]),
                score: 0.0,
            },
            Result {
                x: MVector::from_row_slice(3, &[2.0, 3.5, 2.0]),
                score: 0.0,
            },
        ];

        let cent = centroid(&x);
        let ans = vec![1.0, 2.0, 1.0];

        for i in 0..3 {
            assert!((ans[i] - cent.m[i]).abs() < std::f64::EPSILON);
        }
    }

    #[test]
    fn test_nelder_sin() {
        let f = |x: &[f64]| (x[0].sin() * x[1].cos()) * (1.0 / (x[2].abs() + 1.0));

        let x = vec![
            Result {
                x: MVector::from_row_slice(3, &[0.0f64, 0.0, 0.0]),
                score: 0.0,
            },
            Result {
                x: MVector::from_row_slice(3, &[1.0f64, 1.5, 1.0]),
                score: 0.0,
            },
            Result {
                x: MVector::from_row_slice(3, &[2.0f64, 2.5, 2.0]),
                score: 0.0,
            },
        ];
        info!("centroid: {:?}", centroid(&x));

        let mut x = vec![1.0f64, 2.0, 3.0];

        let f_x = f(&x);
        let result = nelder_mead(
            &f,
            &mut x,
            &Nelder {
                step: 0.1,
                alpha: 1.0,
                gamma: 2.0,
                rho: 0.5,
                sigma: 0.5,
            },
            1e-6,
            100,
        );

        info!("Start result: {}", f_x);
        info!("Optimimum result: {:?} at x: {:?}", result, x);
        assert!((result + 1.0).abs() < 0.01);
    }

    #[test]
    fn test_nelder_rosenbrock() {
        let mut x = vec![2.0f64, 3.0];

        let f_x = rosenbrock2d(&x);
        let result = nelder_mead(
            &rosenbrock2d,
            &mut x,
            &Nelder {
                step: 0.1,
                alpha: 1.0,
                gamma: 2.0,
                rho: 0.5,
                sigma: 0.5,
            },
            1e-6,
            100,
        );

        info!("Start result: {}", f_x);
        info!("Optimimum result: {:?} at x: {:?}", result, x);

        assert!((result).abs() < 0.01);
        assert!(x.iter().fold(0.0, |s, x| s + (x - 1.0).powi(2)) < 0.001);
    }

    fn rosenbrock2d(x: &[f64]) -> f64 {
        (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0].powi(2)).powi(2)
    }

    #[test]
    fn test_nelder_himmelblau() {
        let x_all = vec![[5.0, -5.0], [-5.0, 5.0], [-5.0, -5.0], [0.0, 0.0]];
        let ans_all = vec![
            [3.584_428, -1.848_126],
            [-2.805_118, 3.131_312],
            [-3.779_310, -3.283_186],
            [3.0, 2.0],
        ];

        for (x, ans) in x_all.iter().zip(ans_all) {
            let f_x = himmelblau(x);
            let mut y = *x;
            let result = nelder_mead(
                &himmelblau,
                &mut y,
                &Nelder {
                    step: 0.1,
                    alpha: 1.0,
                    gamma: 2.0,
                    rho: 0.5,
                    sigma: 0.5,
                },
                1e-10,
                100,
            );

            info!("Start result: {}", f_x);
            info!("Optimimum result: {:?} at x: {:?}", result, y);

            assert!((result).abs() < 1e-5);
            let vec_error = y.iter()
                .zip(ans.iter())
                .fold(0.0, |s, (x, a)| s + (x - a).powi(2));
            info!(
                "vecs error {}, answer {:?} obtained {:?}",
                vec_error, ans, y
            );
            assert!(vec_error < 1e-5);
        }
    }

    fn himmelblau(x: &[f64]) -> f64 {
        // this has one local maximum at f(-0.270845, -0.923039) = 181.617
        // and four identical local minima at
        // f(3.0, 2.0) = 0.0
        // f(-2.805118, 3.131312) = 0.0
        // f(-3.779310, -3.283186) = 0.0
        // f(3.584428, -1.848126) = 0.0

        (x[0].powi(2) + x[1] - 11.0).powi(2) + (x[0] + x[1].powi(2) - 7.0).powi(2)
    }

}
