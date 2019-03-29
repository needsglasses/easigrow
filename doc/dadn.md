# Predicting the rate of fatigue crack growth

The following equations are build into the program. To add a new
equation requires adding a new subroutine such as those in the
`dadn.rs` file.

The Paris equation is
$$
\begin{equation}
\frac{da}{dN} = C_{\rm p} \Delta K^{m_{\rm p}}
\label{paris-eqn}
\end{equation}
$$
where $C_{\rm p}$ and $m_{\rm p}$ are constants related to the material, loading and
environment. We use the subscript $p$ to indicate that they are for the
Paris-Erdogen equation. It should be noted that although there are
many similar equations that have $C_{\rm p}$ and $m_{\rm p}$ like coefficients, they cannot
be individually used in other equations without obtaining all the best fit
coefficients again.

