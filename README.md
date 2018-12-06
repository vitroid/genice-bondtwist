# genice-svg

A [GenIce](https://github.com/vitroid/GenIce) plugin to analyze bond twists.

## Requirements

* [GenIce](https://github.com/vitroid/GenIce) >=0.23.

## Installation from PyPI

    % pip install genice-bondtwist

## Manual Installation

### System-wide installation

    % make install

### Private installation

Copy the files in genice_bondtwist/formats/ into your local formats folder of GenIce.

## Usage

	% genice T2 -f bondtwist > T2.btwc

## Test in place

    % make test
