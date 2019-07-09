package psqr

import (
	"sort"
)

type PSQuantile struct {
	s []*Quantile
	m map[float64]int
}

func NewPSQuantile(f ...float64) *PSQuantile {
	q := &PSQuantile{
		m: make(map[float64]int),
	}

	q.Init(f...)

	return q
}

func (s *PSQuantile) Init(f ...float64) {
	for _, v := range f {
		idx := len(s.s)
		s.s = append(s.s, NewQuantile(v))
		s.m[v] = idx
	}
}

func (s *PSQuantile) Append(val float64) {
	for _, v := range s.s {
		v.Append(val)
	}
}

func (s *PSQuantile) Quantile(val float64) float64 {
	return s.s[s.m[val]].Value()
}


// P-square maitains five markers that store points.
const nMarkers = 5

// Quantile represents an estimated p-quantile of a stream of observations.
type Quantile struct {
	p      float64
	filled bool
	cnt    int

	// marker positions, 1..nMarkers
	pos [nMarkers]int
	// desired marker positions
	npos [nMarkers]float64
	// increament in desired marker positions
	dn [nMarkers]float64
	// marker heights that store observations
	heights [nMarkers]float64
}

// NewQuantile returns new p-quantile.
func NewQuantile(p float64) *Quantile {
	if p < 0 || p > 1 {
		panic("p-quantile is out of range")
	}
	q := &Quantile{
		p: p,
	}
	q.Reset()
	return q
}

// Reset resets the quantile.
func (q *Quantile) Reset() {
	p := q.p
	q.filled = false
	q.cnt = 0
	for i := 0; i < len(q.pos); i++ {
		q.pos[i] = i
	}
	q.npos = [...]float64{
		0,
		2 * p,
		4 * p,
		2 + 2*p,
		4,
	}
	q.dn = [...]float64{
		0,
		p / 2,
		p,
		(1 + p) / 2,
		1,
	}
}

// Append appends v to the stream of observations.
func (q *Quantile) Append(v float64) {
	if q.cnt != nMarkers {
		// no required number of observations has been appended yet
		q.heights[q.cnt] = v
		q.cnt++
		return
	}
	if !q.filled {
		q.filled = true
		sort.Float64s(q.heights[:q.cnt])
	}
	q.append(v)

}

func (q *Quantile) append(v float64) {
	l := q.cnt - 1

	k := -1
	if v < q.heights[0] {
		k = 0
		q.heights[0] = v
	} else if q.heights[l] <= v {
		k = l - 1
		q.heights[l] = v
	} else {
		for i := 1; i <= l; i++ {
			if q.heights[i-1] <= v && v < q.heights[i] {
				k = i - 1
				break
			}
		}
	}

	for i := 0; i < nMarkers; i++ {
		// increment positions greater than k
		if i > k {
			q.pos[i]++
		}
		// update desired positions for all markers
		q.npos[i] += q.dn[i]
	}

	q.adjustHeights()
}

func (q *Quantile) adjustHeights() {
	for i := 1; i < len(q.heights)-1; i++ {
		d := int(q.npos[i]) - q.pos[i]

		ni1 := q.pos[i+1] - q.pos[i]
		ni2 := q.pos[i] - q.pos[i-1]

		n1 := float64(ni1)
		n2 := float64(ni2)

		h := q.heights[i]
		hp1 := q.heights[i+1]
		hm1 := q.heights[i-1]

		h1 := hp1 - h
		h2 := h - hm1

		z1 := h1 / n1
		z2 := h2 / n2

		if d >= 1 && ni1 > 1 {
			b1 := (n2 + 1) * z1
			b2 := (n1 - 1) * z2
			hi := h + (b1+b2)/(n1+n2)

			if hm1 < hi && hi < hp1 {
				q.heights[i] = hi
			} else {
				// use linear formula
				q.heights[i] = h + z1
			}

			q.pos[i]++
		} else if d <= -1 && ni2 > -1 {
			b1 := (n2 - 1) * z1
			b2 := (n1 + 1) * z2
			hi := h - (b1+b2)/(n1+n2)

			if hm1 < hi && hi < hp1 {
				q.heights[i] = hi
			} else {
				// use linear formula
				q.heights[i] = h - z2
			}

			q.pos[i]--
		}
	}
}

// Value returns the current estimate of p-quantile.
func (q *Quantile) Value() float64 {
	if !q.filled {
		// a fast path when not enought observations has been stored yet
		switch q.cnt {
		case 0:
			return 0
		case 1:
			return q.heights[0]
		}
		sort.Float64s(q.heights[:q.cnt])
		rank := int(q.p * float64(q.cnt))
		return q.heights[rank]
	}
	// if initialised with nMarkers observations third height stores current
	// estimate of p-quantile
	return q.heights[2]
}
