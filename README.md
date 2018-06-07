# genice-svg

A [GenIce](https://github.com/vitroid/GenIce) plugin to illustrate the structure in SVG format.

## Requirements

* [GenIce](https://github.com/vitroid/GenIce) >=0.23.
* svgwrite.

## Installation

### System-wide installation

    % make install

### Private installation

Copy the files in genice_svg/formats/ into your local formats folder of GenIce.

## Usage

	% genice CS2 -r 3 3 3 -f svg_poly > CS2.svg

## Test in place

    % make test
