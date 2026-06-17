# plot_corner_pdf.m
This MATLAB function creates a corner plot of the discretized PDF represented by $X$ and $P$, where $X$ is the state coordinate matrix and $P$ is the probability at each of these coordinate points. <br>

$$
\begin{gather}
    X \in \mathbb{R}^{n\times d} \\
    P \in \mathbb{R}^{n\times 1} 
\end{gather}
$$

This function is meant to be used in tandem with [plot_nongaussian_surface.m](https://github.com/bhanson10/plot_nongaussian_surface) and [plot_gaussian_ellipsoid.m](https://github.com/bhanson10/plot_gaussian_ellipsoid). <br>

Please direct any questions to blhanson@ucsd.edu. <br><br>

![plot_corner_pdf](plot_corner_pdf.png)
