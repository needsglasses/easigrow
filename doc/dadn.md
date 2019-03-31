# Predicting the rate of fatigue crack growth

The following equations are build into the program. To add a new
equation requires adding a new subroutine such as those in the
`dadn.rs` file.

The Paris equation is

$\frac{da}{dN} = C_p \Delta K^{m_p}$

where $C_p$ and $m_p$ are constants related to the material, loading and
environment. We use the subscript $p$ to indicate that they are for the
Paris-Erdogen equation. It should be noted that although there are
many similar equations that have $C_p$ and $m_p$ like coefficients, they cannot
be individually used in other equations without obtaining all the best fit
coefficients again.

