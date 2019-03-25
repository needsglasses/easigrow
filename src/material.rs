//! Database of fatigue parameters for typical materials.
use dadn;

#[derive(Debug, Clone)]
pub struct Properties {
    pub yield_stress: f64,
    pub k1c: f64,
    pub youngs_modulus: f64,
}

pub struct EqnProperty {
    pub prop: Props,
    pub name: &'static str,
    pub cite: &'static str,
    pub units: &'static str,
    pub eqn: Box<dadn::DaDn>,
}

#[derive(Debug, Clone)]
pub enum Props {
    YieldStress,
    YoungsModulus,
    PoissonRatio,
    Kc,
    K1c,
    Dadn,
}

#[derive(Debug, Clone)]
pub struct Property {
    pub prop: Props,
    pub cite: &'static str,
    pub value: f64,
    pub units: &'static str,
}

pub fn get_all_dadns() -> [EqnProperty; 33] {
    let aluminium = [
        Property {
            prop: Props::YieldStress,
            cite: "",
            value: 450.0,
            units: "MPa",
        },
        Property {
            prop: Props::YoungsModulus,
            cite: "",
            value: 71000.0,
            units: "MPa",
        },
        Property {
            prop: Props::K1c,
            cite: "",
            value: 31.0,
            units: "MPa sqrt(m)",
        },
    ];

    let _x = [aluminium];
    let _cite = "".to_string();
    [
        EqnProperty {
            prop: Props::Dadn,
            name: "paris:default",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Paris::new(&vec![1.00e-10, 3.0], "".to_string()))
                as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "walker:default",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Walker::new(&vec![1.00e-10, 0.5, 3.0], "".to_string()))
                as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:default",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![1.00e-10, 3.0, 60.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "nasgro:default",
            cite: "[nasgro:aa7050-t7451-LT, NASGR04.0]",
            units: "m",
            eqn: Box::new(dadn::Nasgro::new(
                &vec![
                    0.3, 2.0, 35.16, 0.80, 2.20, 0.1, 1.0, 1.0, 6.35e-10, 2.50, 38.1e-6
                ],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "nasgro:aa7050-t7451-LT",
            cite: "[NASGR04.0]",
            // // coeffs: vec![0.3, 2.0, 32.0, 0.80, 2.20, 0.1, 1.0, 1.0, 0.25e-7, 2.50, 1.5e-3], // in, ksi.sqrt(in)
            units: "m",
            eqn: Box::new(dadn::Nasgro::new(
                &vec![
                    0.3, 2.0, 35.16, 0.80, 2.20, 0.1, 1.0, 1.0, 6.35e-10, 2.50, 38.1e-6
                ],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "burchill:default",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Burchill::new(
                &vec![1.00e-10, 3.0, 1.00e-10, 3.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "hartman:default",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Hartman::new(
                &vec![1.00e-10, 1.0, 30.0, 3.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "paris:newman-7050t7451",
            cite: "[none]",
            units: "m",
            eqn: Box::new(dadn::Paris::new(&vec![1.593e-11, 3.668], "".to_string()))
                as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t3-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![7.13e-9, 2.70, 71.3],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t351-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![5.00e-9, 2.88, 63.2],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t4-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![8.57e-9, 2.60, 58.1],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t6-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.00e-8, 2.62, 69.8],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t8-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![1.33e-8, 2.65, 65.3],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2024t851-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![7.72e-9, 2.78, 61.4],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2219t851-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![4.84e-8, 2.16, 57.5],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa2618t6-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![8.56e-9, 2.58, 45.9],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa6061t6-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.27e-7, 60.1, 1.66],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa6061t651-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![9.60e-8, 1.84, 41.2],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7010t73651-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.06e-8, 2.46, 46.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7050t7352-forging",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.75e-9, 3.29, 64.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7050t73651-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![4.11e-9, 2.98, 55.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7075t6-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![1.37e-8, 3.02, 63.9],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7075t7351",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![6.27e-9, 2.78, 55.8],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7175t3652-forging",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.61e-9, 2.91, 38.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7178t651-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![3.74e-8, 2.06, 30.7],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7475t7351-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![3.24e-8, 2.32, 78.2],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7475t76-sheet",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![6.54e-8, 2.18, 79.9],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:aa7475t7651-plate",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![9.30e-9, 2.73, 63.1],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:a357t6-sandcasting",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![2.19e-9, 2.94, 41.5],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "forman:a357t6-investmentcasting",
            cite: "[Schwarmann86]",
            units: "m",
            eqn: Box::new(dadn::Forman::new(
                &vec![6.65e-9, 2.40, 38.2],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "hartman:jones13-aa7050t7451",
            cite: "[jones13]",
            units: "m",
            eqn: Box::new(dadn::Hartman::new(
                &vec![7.0e-10, 0.1, 47.0, 2.0],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "white:barter14-aa7050t7451",
            cite: "[white15]",

            units: "m",
            eqn: Box::new(dadn::White::new(
                &vec![
                    0.25481858f64,
                    1.10247048,
                    4.35831677,
                    23.08586582,
                    0.03420171,
                    0.4717843,
                    31.54,
                ],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
        EqnProperty {
            prop: Props::Dadn,
            name: "white:chan16-aa7050t7451",
            cite: "[]",

            units: "m",
            eqn: Box::new(dadn::White::new(
                &vec![
                    0.2918621458996122,
                    1.2635107616404941,
                    3.5528305197144334,
                    22.243180576185246,
                    0.03924086425080324,
                    0.5551311271413691,
                    41.45917108627669,
                ],
                "".to_string(),
            )) as Box<dadn::DaDn>,
        },
    ]
}

#[cfg(test)]
mod tests {
    use table;

    #[test]
    fn check_pairtable() {
        // ID Code = M7GJ11AC1
        let nasgro_table = table::PairTable::new(
            vec![0.08, 0.1, 0.4, 0.5, 0.7, 0.8],
            vec![
                vec![
                    6.8865, 7.1373, 8.5740, 11.2303, 12.6427, 13.7598, 18.5817, 21.4952, 23.1962
                ],
                vec![
                    3.9084, 3.9241, 4.2994, 4.8354, 5.5089, 6.4799, 7.5038, 8.9768, 10.6903,
                    13.0444, 15.4584, 18.3955, 21.1440, 23.9296, 26.0453, 28.3481, 29.8538,
                ],
                vec![
                    2.5410, 2.5675, 3.0591, 3.6212, 4.1691, 4.7674, 5.3145, 6.1336, 7.2176, 8.8385,
                    10.6337, 13.1499, 15.5498, 16.7494,
                ],
                vec![
                    4.7534, 4.8525, 5.4633, 6.2946, 7.6432, 9.2110, 11.4217, 13.7010, 16.4040,
                    18.8969, 22.2837, 26.4682, 26.4850,
                ],
                vec![
                    1.9011, 1.9119, 2.2139, 2.8924, 3.3723, 3.8130, 4.2446, 4.8279, 5.4787, 6.1235
                ],
                vec![
                    1.7701, 2.0983, 2.5583, 2.9669, 3.3785, 3.8713, 4.7429, 5.5754, 5.8345
                ],
            ],
            vec![
                vec![
                    8.1560e-07, 1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05,
                    1.8000e-05, 3.2000e-05, 6.0000e-05,
                ],
                vec![
                    5.8598e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05, 1.8000e-05,
                    3.2000e-05, 6.0000e-05, 1.0000e-04, 1.8000e-04, 2.6886e-04,
                ],
                vec![
                    5.8374e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05, 1.8000e-05,
                    3.2000e-05, 4.3702e-05,
                ],
                vec![
                    8.9219e-07, 1.0000e-06, 1.8000e-06, 3.2000e-06, 6.0000e-06, 1.0000e-05,
                    1.8000e-05, 3.2000e-05, 6.0000e-05, 1.0000e-04, 1.8000e-04, 3.2000e-04,
                    3.2064e-04,
                ],
                vec![
                    5.8414e-08, 6.0000e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07,
                    1.0000e-06, 1.8000e-06, 3.2000e-06, 5.4116e-06,
                ],
                vec![
                    6.3591e-08, 1.0000e-07, 1.8000e-07, 3.2000e-07, 6.0000e-07, 1.0000e-06,
                    1.8000e-06, 3.2000e-06, 3.9566e-06,
                ],
            ],
        );

        assert_eq!(nasgro_table.interp(6.8865, 0.08), 8.156e-7);
        assert_eq!(nasgro_table.interp(5.8345, 0.8), 3.9566e-6);
        assert!((nasgro_table.interp(23.1962, 0.08) - 6.0e-5).abs() < 1.0e-7);
        assert!((nasgro_table.interp(4.7534, 0.5) - 8.9219e-7).abs() < 1.0e-7);
        assert!((nasgro_table.interp(26.4850, 0.5) - 3.2064e-4).abs() < 1.0e-7);
    }

}
