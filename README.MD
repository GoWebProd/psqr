# psqr
Go implementation of [The P-Square Algorithm for Dynamic Calculation of Quantiles and Histograms Without Storing Observations][1].

> The algorithm is proposed for dynamic calculation of [..] quantiles. The estimates are produced dynamically as the observations are generated. The observations are not stored, therefore, the algorithm has a very small and fixed storage requirement regardless of the number of observations.

[1]: http://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf

## Usage example

```go
package main

import (
	"fmt"
	
	"github.com/GoWebProd/psqr"
)

func main() {
	pp := NewPSQuantile(0.5, 0.9, 0.99, 0.999)

	pp.Append(0.013163238)
	pp.Append(0.711542201)
	pp.Append(-2.131796046)
	pp.Append(0.244640008)
	pp.Append(-0.211374733)
	pp.Append(4.493872061)

	fmt.Println("0.5 =", pp.Quantile(0.5))
	fmt.Println("0.9 =", pp.Quantile(0.9))
	fmt.Println("0.99 =", pp.Quantile(0.99))
	fmt.Println("0.999 =", pp.Quantile(0.999))
}
```
