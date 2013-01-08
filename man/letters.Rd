\name{letters}
\docType{data}
\alias{letters}
\alias{lets}

\title{Letter Recognition Data Set}

\description{
  The objective is to identify each of a large number of black-and-white rectangular pixel displays as one of the 26 capital letters in the English alphabet. 
}
\usage{data(letters)}

\format{A data frame of 17 variables and 20000 entries.
Description of the variables: 
1. lettr capital letter (26 values from A to Z)
2. x-box horizontal position of box (integer)
3. y-box vertical position of box (integer)
4. width width of box (integer)
5. high height of box (integer)
6. onpix total # on pixels (integer)
7. x-bar mean x of on pixels in box (integer)
8. y-bar mean y of on pixels in box (integer)
9. x2bar mean x variance (integer)
10. y2bar mean y variance (integer)
11. xybar mean x y correlation (integer)
12. x2ybr mean of x * x * y (integer)
13. xy2br mean of x * y * y (integer)
14. x-ege mean edge count left to right (integer)
15. xegvy correlation of x-ege with y (integer)
16. y-ege mean edge count bottom to top (integer)
17. yegvx correlation of y-ege with x (integer)
}

\source{UCI Machine Learning Repository, Creator: David J. Slate, Odesta Corporation; 1890 Maple Ave; Suite 115; Evanston, IL 60201}

\keyword{datasets}
