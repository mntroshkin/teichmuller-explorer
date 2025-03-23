# teichmuller-explorer
A small Python script for generating images of hyperbolic tilings.

### Mathematical Background

It is well-known that a sphere with $g$ handles is glued from a $4g$-gon by pairwise gluing of the sides. Let us take (for $g \geq 2$) a $4g$-gon in the hyperbolic plane, whose identified pairs of sides are equal in length, and the sum of the angles is $2\pi$. Then the resulting surface has a hyperbolic metric; conversely, any hyperbolic surface can be cut (along geodesic segments) into such a $4g$-gon, whose vertices are be glued together into any point chosen on the surface.

If we now consider the universal covering of the surface by the hyperbolic plane, then the corresponding $4g$-gon will tile the plane with its copies. Thus, we can depict hyperbolic metrics on the surface (points in the Teichmüller space) using hyperbolic tilings.

This is a program that draws a fragment of such a random tiling for $g=2$, i.e. using octagons.

### Prerequisites

* Python 3.x
* Python libraries `numpy`, `scipy`, `drawsvg`

### Installation

1. **Clone the repository:**
   
   ```bash
   git clone https://github.com/mntroshkin/teichmuller-explorer.git
   cd my-python-program
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

### Usage

To run the program, execute the file `main.py`. It will output the SVG image of a random hyperbolic tiling by octagons to a file `tiling.svg` in the current folder.

### Acknowledgements

I relied heavily on archived materials from a similar project, which had a greater interactive functional: "The Teichmüller Navigator", created in 1993 by Deva Van Der Werf, David Ben-Zvi, and Paul Burchard. Miraculously the project page has survived to this day: http://www.geom.uiuc.edu/apps/teich-nav/report/report.html.
