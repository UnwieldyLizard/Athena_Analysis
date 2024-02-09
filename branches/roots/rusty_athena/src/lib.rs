use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayView4};
use numpy::{PyArray1, PyArray2, PyArray4, IntoPyArray};
use pyo3::prelude::*;

#[pyfunction]
fn volume_integral(
    q: &PyArray4<f64>,
    zeroth_block_coordinates: &PyArray2<f64>,
    first_block_coordinates: &PyArray2<f64>,
    second_block_coordinates: &PyArray2<f64>,
    zeroth_block_differentials: &PyArray4<f64>,
    first_block_differentials: &PyArray4<f64>,
    second_block_differentials: &PyArray4<f64>,
    zeroth_grid_coordinates: &PyArray1<f64>,
    first_grid_coordinates: &PyArray1<f64>,
    second_grid_coordinates: &PyArray1<f64>,
    zeroth_grid_differentials: &PyArray1<f64>,
    first_grid_differentials: &PyArray1<f64>,
    second_grid_differentials: &PyArray1<f64>,
    integration_order_str: [&str; 3],
    bounds: ((f64, f64), (f64, f64), (f64, f64)),
) -> PyResult<f64> {
    let mut integration_order: [Coordinate; 3] = [Coordinate::X; 3];

    for idx in 0..integration_order_str.len() {
        integration_order[idx] = match integration_order_str[idx] {
            "x" | "X" => Coordinate::X,
            "y" | "Y" => Coordinate::Y,
            "z" | "Z" => Coordinate::Z,
            "r" | "R" => Coordinate::R,
            "phi" | "Phi" => Coordinate::Phi,
            "theta" | "Theta" => Coordinate::Theta,
            _ => panic!("Unknown variables in coordinate system"),
        };
    }

    // Get the ndarray as a Rust array view
    let q_binding = q.readonly();
    let q: ArrayView4<f64> = q_binding.as_array();

    let zeroth_block_coord_binds = zeroth_block_coordinates.readonly();
    let zeroth_block_coordinates: ArrayView2<f64> = zeroth_block_coord_binds.as_array();
    let first_block_coord_binds = first_block_coordinates.readonly();
    let first_block_coordinates: ArrayView2<f64> = first_block_coord_binds.as_array();
    let second_block_coord_binds = second_block_coordinates.readonly();
    let second_block_coordinates: ArrayView2<f64> = second_block_coord_binds.as_array();

    let zeroth_grid_coord_binds = zeroth_grid_coordinates.readonly();
    let zeroth_grid_coordinates: ArrayView1<f64> = zeroth_grid_coord_binds.as_array();
    let first_grid_coord_binds = first_grid_coordinates.readonly();
    let first_grid_coordinates: ArrayView1<f64> = first_grid_coord_binds.as_array();
    let second_grid_coord_binds = second_grid_coordinates.readonly();
    let second_grid_coordinates: ArrayView1<f64> = second_grid_coord_binds.as_array();

    let zeroth_grid_diff_binds = zeroth_grid_differentials.readonly();
    let zeroth_grid_differentials: ArrayView1<f64> = zeroth_grid_diff_binds.as_array();
    let first_grid_diff_binds = first_grid_differentials.readonly();
    let first_grid_differentials: ArrayView1<f64> = first_grid_diff_binds.as_array();
    let second_grid_diff_binds = second_grid_differentials.readonly();
    let second_grid_differentials: ArrayView1<f64> = second_grid_diff_binds.as_array();

    let zeroth_block_diff_binds = zeroth_block_differentials.readonly();
    let first_block_diff_binds = first_block_differentials.readonly();
    let second_block_diff_binds = second_block_differentials.readonly();
    let first_block_diff: ArrayView4<f64> = match integration_order[0] {
        Coordinate::Z | Coordinate::Theta => zeroth_block_diff_binds.as_array(),
        Coordinate::Y | Coordinate::Phi => first_block_diff_binds.as_array(),
        Coordinate::X | Coordinate::R => second_block_diff_binds.as_array(),
    };

    let first_int_q = initial_integral(
        &q,
        &zeroth_block_coordinates,
        &first_block_coordinates,
        &second_block_coordinates,
        &first_block_diff,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        bounds.0,
        &integration_order[0],
    );

    let second_int_q = second_integral(
        &first_int_q,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        &zeroth_grid_differentials,
        &first_grid_differentials,
        &second_grid_differentials,
        bounds.1,
        &integration_order,
    );

    let final_int_q = final_integral(
        &second_int_q,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        &zeroth_grid_differentials,
        &first_grid_differentials,
        &second_grid_differentials,
        bounds.2,
        &integration_order,
    );

    // Convert the result back to a PyArrayDyn
    //final_int_q.into_pyarray(py)
    Ok(final_int_q)
}

#[pyfunction]
fn area_integral<'a>(
    py: Python<'a>,
    q: &PyArray4<f64>,
    zeroth_block_coordinates: &PyArray2<f64>,
    first_block_coordinates: &PyArray2<f64>,
    second_block_coordinates: &PyArray2<f64>,
    zeroth_block_differentials: &PyArray4<f64>,
    first_block_differentials: &PyArray4<f64>,
    second_block_differentials: &PyArray4<f64>,
    zeroth_grid_coordinates: &PyArray1<f64>,
    first_grid_coordinates: &PyArray1<f64>,
    second_grid_coordinates: &PyArray1<f64>,
    zeroth_grid_differentials: &PyArray1<f64>,
    first_grid_differentials: &PyArray1<f64>,
    second_grid_differentials: &PyArray1<f64>,
    integration_order_str: [&str; 2],
    bounds: ((f64, f64), (f64, f64)),
) -> PyResult<&'a PyArray1<f64>> {
    let mut integration_order: [Coordinate; 3] = [Coordinate::X; 3];

    for idx in 0..integration_order_str.len() {
        integration_order[idx] = match integration_order_str[idx] {
            "x" | "X" => Coordinate::X,
            "y" | "Y" => Coordinate::Y,
            "z" | "Z" => Coordinate::Z,
            "r" | "R" => Coordinate::R,
            "phi" | "Phi" => Coordinate::Phi,
            "theta" | "Theta" => Coordinate::Theta,
            _ => panic!("Unknown variables in coordinate system"),
        };
    }

    // Get the ndarray as a Rust array view
    let q_binding = q.readonly();
    let q: ArrayView4<f64> = q_binding.as_array();

    let zeroth_block_coord_binds = zeroth_block_coordinates.readonly();
    let zeroth_block_coordinates: ArrayView2<f64> = zeroth_block_coord_binds.as_array();
    let first_block_coord_binds = first_block_coordinates.readonly();
    let first_block_coordinates: ArrayView2<f64> = first_block_coord_binds.as_array();
    let second_block_coord_binds = second_block_coordinates.readonly();
    let second_block_coordinates: ArrayView2<f64> = second_block_coord_binds.as_array();

    let zeroth_grid_coord_binds = zeroth_grid_coordinates.readonly();
    let zeroth_grid_coordinates: ArrayView1<f64> = zeroth_grid_coord_binds.as_array();
    let first_grid_coord_binds = first_grid_coordinates.readonly();
    let first_grid_coordinates: ArrayView1<f64> = first_grid_coord_binds.as_array();
    let second_grid_coord_binds = second_grid_coordinates.readonly();
    let second_grid_coordinates: ArrayView1<f64> = second_grid_coord_binds.as_array();

    let zeroth_grid_diff_binds = zeroth_grid_differentials.readonly();
    let zeroth_grid_differentials: ArrayView1<f64> = zeroth_grid_diff_binds.as_array();
    let first_grid_diff_binds = first_grid_differentials.readonly();
    let first_grid_differentials: ArrayView1<f64> = first_grid_diff_binds.as_array();
    let second_grid_diff_binds = second_grid_differentials.readonly();
    let second_grid_differentials: ArrayView1<f64> = second_grid_diff_binds.as_array();

    let zeroth_block_diff_binds = zeroth_block_differentials.readonly();
    let first_block_diff_binds = first_block_differentials.readonly();
    let second_block_diff_binds = second_block_differentials.readonly();
    let first_block_diff: ArrayView4<f64> = match integration_order[0] {
        Coordinate::Z | Coordinate::Theta => zeroth_block_diff_binds.as_array(),
        Coordinate::Y | Coordinate::Phi => first_block_diff_binds.as_array(),
        Coordinate::X | Coordinate::R => second_block_diff_binds.as_array(),
    };

    let first_int_q = initial_integral(
        &q,
        &zeroth_block_coordinates,
        &first_block_coordinates,
        &second_block_coordinates,
        &first_block_diff,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        bounds.0,
        &integration_order[0],
    );

    let second_int_q = second_integral(
        &first_int_q,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        &zeroth_grid_differentials,
        &first_grid_differentials,
        &second_grid_differentials,
        bounds.1,
        &integration_order,
    );

    // Convert the result back to a PyArrayDyn
    let final_int_q = second_int_q.into_pyarray(py);
    Ok(final_int_q)
}

#[pyfunction]
fn single_integral<'a>(
    py: Python<'a>,
    q: &PyArray4<f64>,
    zeroth_block_coordinates: &PyArray2<f64>,
    first_block_coordinates: &PyArray2<f64>,
    second_block_coordinates: &PyArray2<f64>,
    zeroth_block_differentials: &PyArray4<f64>,
    first_block_differentials: &PyArray4<f64>,
    second_block_differentials: &PyArray4<f64>,
    zeroth_grid_coordinates: &PyArray1<f64>,
    first_grid_coordinates: &PyArray1<f64>,
    second_grid_coordinates: &PyArray1<f64>,
    integration_order_str: [&str; 1],
    bounds: (f64, f64),
) -> PyResult<&'a PyArray2<f64>> {
    let mut integration_order: [Coordinate; 3] = [Coordinate::X; 3];

    for idx in 0..integration_order_str.len() {
        integration_order[idx] = match integration_order_str[idx] {
            "x" | "X" => Coordinate::X,
            "y" | "Y" => Coordinate::Y,
            "z" | "Z" => Coordinate::Z,
            "r" | "R" => Coordinate::R,
            "phi" | "Phi" => Coordinate::Phi,
            "theta" | "Theta" => Coordinate::Theta,
            _ => panic!("Unknown variables in coordinate system"),
        };
    }

    // Get the ndarray as a Rust array view
    let q_binding = q.readonly();
    let q: ArrayView4<f64> = q_binding.as_array();

    let zeroth_block_coord_binds = zeroth_block_coordinates.readonly();
    let zeroth_block_coordinates: ArrayView2<f64> = zeroth_block_coord_binds.as_array();
    let first_block_coord_binds = first_block_coordinates.readonly();
    let first_block_coordinates: ArrayView2<f64> = first_block_coord_binds.as_array();
    let second_block_coord_binds = second_block_coordinates.readonly();
    let second_block_coordinates: ArrayView2<f64> = second_block_coord_binds.as_array();

    let zeroth_grid_coord_binds = zeroth_grid_coordinates.readonly();
    let zeroth_grid_coordinates: ArrayView1<f64> = zeroth_grid_coord_binds.as_array();
    let first_grid_coord_binds = first_grid_coordinates.readonly();
    let first_grid_coordinates: ArrayView1<f64> = first_grid_coord_binds.as_array();
    let second_grid_coord_binds = second_grid_coordinates.readonly();
    let second_grid_coordinates: ArrayView1<f64> = second_grid_coord_binds.as_array();

    let zeroth_block_diff_binds = zeroth_block_differentials.readonly();
    let first_block_diff_binds = first_block_differentials.readonly();
    let second_block_diff_binds = second_block_differentials.readonly();
    let first_block_diff: ArrayView4<f64> = match integration_order[0] {
        Coordinate::Z | Coordinate::Theta => zeroth_block_diff_binds.as_array(),
        Coordinate::Y | Coordinate::Phi => first_block_diff_binds.as_array(),
        Coordinate::X | Coordinate::R => second_block_diff_binds.as_array(),
    };

    let first_int_q = initial_integral(
        &q,
        &zeroth_block_coordinates,
        &first_block_coordinates,
        &second_block_coordinates,
        &first_block_diff,
        &zeroth_grid_coordinates,
        &first_grid_coordinates,
        &second_grid_coordinates,
        bounds,
        &integration_order[0],
    );

    // Convert the result back to a PyArrayDyn
    let final_int_q = first_int_q.into_pyarray(py);
    Ok(final_int_q)
}

fn initial_integral<'a>(
    q: &ArrayView4<'a, f64>,
    zeroth_block_coordinates: &ArrayView2<'a, f64>,
    first_block_coordinates: &ArrayView2<'a, f64>,
    second_block_coordinates: &ArrayView2<'a, f64>,
    block_differential: &ArrayView4<'a, f64>,
    zeroth_grid_coordinates: &ArrayView1<'a, f64>,
    first_grid_coordinates: &ArrayView1<'a, f64>,
    second_grid_coordinates: &ArrayView1<'a, f64>,
    bounds: (f64, f64),
    coordinate: &Coordinate,
) -> Array2<f64> {
    let remaining_coords: (usize, usize) = match coordinate {
        Coordinate::Z | Coordinate::Theta => (1, 2),
        Coordinate::Y | Coordinate::Phi => (0, 2),
        Coordinate::X | Coordinate::R => (0, 1),
    };

    let remaining_grid_coords = match coordinate {
        Coordinate::Z | Coordinate::Theta => (first_grid_coordinates, second_grid_coordinates),
        Coordinate::Y | Coordinate::Phi => (zeroth_grid_coordinates, second_grid_coordinates),
        Coordinate::X | Coordinate::R => (zeroth_grid_coordinates, first_grid_coordinates),
    };

    let remaining_block_coords = match coordinate {
        Coordinate::Z | Coordinate::Theta => (first_block_coordinates, second_block_coordinates),
        Coordinate::Y | Coordinate::Phi => (zeroth_block_coordinates, second_block_coordinates),
        Coordinate::X | Coordinate::R => (zeroth_block_coordinates, first_block_coordinates),
    };

    let int_block_coords = match coordinate {
        Coordinate::Z | Coordinate::Theta => zeroth_block_coordinates,
        Coordinate::Y | Coordinate::Phi => first_block_coordinates,
        Coordinate::X | Coordinate::R => second_block_coordinates,
    };

    let mut int_q_dx: Array2<f64> =
        Array2::zeros((remaining_grid_coords.0.len(), remaining_grid_coords.1.len()));

    for mesh_idx in 0..q.shape()[0] {
        for idx0 in 0..q.shape()[1] {
            for idx1 in 0..q.shape()[2] {
                for idx2 in 0..q.shape()[3] {
                    let idx_arr: [usize; 3] = [idx0, idx1, idx2];

                    if int_block_coords[[mesh_idx, idx_arr[coordinate.idx()]]] < bounds.0
                        || int_block_coords[[mesh_idx, idx_arr[coordinate.idx()]]] > bounds.1
                    {
                        continue;
                    }
                    
                    let first_coord_idx = binary_search(
                        &remaining_grid_coords.0,
                        &remaining_block_coords.0[[mesh_idx, idx_arr[remaining_coords.0]]],
                    );

                    let second_coord_idx = binary_search(
                        &remaining_grid_coords.1,
                        &remaining_block_coords.1[[mesh_idx, idx_arr[remaining_coords.1]]],
                    );

                    match first_coord_idx {
                        Some(first_coord_idx) => match second_coord_idx {
                            Some(second_coord_idx) => {
                                int_q_dx[[first_coord_idx, second_coord_idx]] += match coordinate {
                                    Coordinate::Phi => {
                                        q[[mesh_idx, idx0, idx1, idx2]]
                                            * block_differential
                                                [[mesh_idx, idx0, idx1, idx2]]
                                            * second_block_coordinates
                                                [[mesh_idx, idx_arr[Coordinate::R.idx()]]]
                                    } // the r in rdphi
                                    Coordinate::Theta => {
                                        q[[mesh_idx, idx0, idx1, idx2]]
                                            * block_differential
                                                [[mesh_idx, idx0, idx1, idx2]]
                                            * second_block_coordinates
                                                [[mesh_idx, idx_arr[Coordinate::R.idx()]]]
                                            * first_block_coordinates
                                                [[mesh_idx, idx_arr[Coordinate::Phi.idx()]]]
                                            .sin()
                                    } // the rsin(phi) in rsin(phi)dtheta
                                    _ => {
                                        q[[mesh_idx, idx0, idx1, idx2]]
                                            * block_differential
                                                [[mesh_idx, idx0, idx1, idx2]]
                                    }
                                }
                            }
                            None => panic!("Value finding error (Second Index)"),
                        },
                        None => panic!("Value finding error (First Index)"),
                    }
                }
            }
        }
    }

    int_q_dx
}

fn second_integral<'a>(
    q: &Array2<f64>,
    zeroth_grid_coordinates: &ArrayView1<'a, f64>,
    first_grid_coordinates: &ArrayView1<'a, f64>,
    second_grid_coordinates: &ArrayView1<'a, f64>,
    zeroth_grid_differentials: &ArrayView1<'a, f64>,
    first_grid_differentials: &ArrayView1<'a, f64>,
    second_grid_differentials: &ArrayView1<'a, f64>,
    bounds: (f64, f64),
    integration_order: &[Coordinate; 3],
) -> Array1<f64> {
    let grid_coord = match integration_order[1] {
        Coordinate::Z | Coordinate::Theta => zeroth_grid_coordinates,
        Coordinate::Y | Coordinate::Phi => first_grid_coordinates,
        Coordinate::X | Coordinate::R => second_grid_coordinates,
    };

    let grid_diff = match integration_order[1] {
        Coordinate::Z | Coordinate::Theta => zeroth_grid_differentials,
        Coordinate::Y | Coordinate::Phi => first_grid_differentials,
        Coordinate::X | Coordinate::R => second_grid_differentials,
    };

    let remaining_grid_coord = match integration_order[2] {
        Coordinate::Z | Coordinate::Theta => zeroth_grid_coordinates,
        Coordinate::Y | Coordinate::Phi => first_grid_coordinates,
        Coordinate::X | Coordinate::R => second_grid_coordinates,
    };

    let mut int_q_dxdy: Array1<f64> = Array1::zeros(remaining_grid_coord.len());

    for idx0 in 0..q.shape()[0] {
        for idx1 in 0..q.shape()[1] {
            let idx_arr: [usize; 3] = match integration_order[0] {
                Coordinate::Z | Coordinate::Theta => [0, idx0, idx1],
                Coordinate::Y | Coordinate::Phi => [idx0, 0, idx1],
                Coordinate::X | Coordinate::R => [idx0, idx1, 0],
            };

            if grid_coord[idx_arr[integration_order[1].idx()]] <= bounds.0 as f64
                || grid_coord[idx_arr[integration_order[1].idx()]] >= bounds.1 as f64
            {
                continue;
            }

            int_q_dxdy[idx_arr[integration_order[2].idx()]] += match integration_order[1] {
                Coordinate::Phi => {
                    q[[idx0, idx1]]
                        * grid_diff[idx_arr[integration_order[1].idx()]]
                        * second_grid_coordinates[idx_arr[Coordinate::R.idx()]]
                }
                _ => q[[idx0, idx1]] * grid_diff[idx_arr[integration_order[1].idx()]],
            }
        }
    }

    int_q_dxdy
}

fn final_integral<'a>(
    q: &Array1<f64>,
    zeroth_grid_coordinates: &ArrayView1<'a, f64>,
    first_grid_coordinates: &ArrayView1<'a, f64>,
    second_grid_coordinates: &ArrayView1<'a, f64>,
    zeroth_grid_differentials: &ArrayView1<'a, f64>,
    first_grid_differentials: &ArrayView1<'a, f64>,
    second_grid_differentials: &ArrayView1<'a, f64>,
    bounds: (f64, f64),
    integration_order: &[Coordinate; 3],
) -> f64 {
    let grid_coord = match integration_order[2] {
        Coordinate::Z | Coordinate::Theta => zeroth_grid_coordinates,
        Coordinate::Y | Coordinate::Phi => first_grid_coordinates,
        Coordinate::X | Coordinate::R => second_grid_coordinates,
    };

    let grid_diff = match integration_order[2] {
        Coordinate::Z | Coordinate::Theta => zeroth_grid_differentials,
        Coordinate::Y | Coordinate::Phi => first_grid_differentials,
        Coordinate::X | Coordinate::R => second_grid_differentials,
    };

    let mut int_q_dxdydz: f64 = 0.0;

    for idx0 in 0..q.len() {
        if grid_coord[idx0] <= bounds.0 as f64 || grid_coord[idx0] >= bounds.1 as f64 {
            continue;
        }

        int_q_dxdydz += q[idx0] * grid_diff[idx0];
    }

    int_q_dxdydz
}

#[derive(Debug, Copy, Clone)]
enum Coordinate {
    X,
    Y,
    Z,
    R,
    Phi,
    Theta, //Uses math convention so phi differential is same in cyl and sph coords
}

impl Coordinate {
    fn idx(&self) -> usize {
        match self {
            Coordinate::X | Coordinate::R => 2,
            Coordinate::Y | Coordinate::Phi => 1,
            Coordinate::Z | Coordinate::Theta => 0,
        }
    }

    fn idx_tuple<'a, T>(&self, tuple: (&'a T, &'a T, &'a T)) -> &'a T {
        match self {
            Coordinate::X | Coordinate::R => tuple.2,
            Coordinate::Y | Coordinate::Phi => tuple.1,
            Coordinate::Z | Coordinate::Theta => tuple.0,    
        }
    }
}

fn binary_search(a: &ArrayView1<f64>, target_value: &f64) -> Option<usize> {
    let mut low: i16 = 0;
    let mut high: i16 = a.len() as i16 - 1;

    while low <= high {
        let mid = ((high - low) / 2) + low;
        let mid_index = mid as usize;
        let val = &a[mid_index];
        
        if val == target_value {
            return Some(mid_index);
        }

        // Search values that are greater than val - to right of current mid_index
        if val < target_value {
            low = mid + 1;
        }

        // Search values that are less than val - to the left of current mid_index
        if val > target_value {
            high = mid - 1;
        }
    }
    None
}

/// A Python module implemented in Rust.
#[pymodule]
fn rusty_athena(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(volume_integral, m)?)?;
    m.add_function(wrap_pyfunction!(area_integral, m)?)?;
    m.add_function(wrap_pyfunction!(single_integral, m)?)?;
    Ok(())
}
