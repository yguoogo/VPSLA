# VPSLA
## Authors: Hayden Brochu, Guanxu Yu


VPSLA is a visualization and prioritization tool of suboptimal local alignments

##  Package features
- Render an efficient way to find suboptimal local alignments
- Allow user to set the threshold.
- Enable customized affine indel penalty
- Use BLOSUM62 matrix to calculate mismatch penalty
- Repot all suboptimal local alignments result
- Plot alignment graph, allow users to view the entire set of alignments.

## Prerequisite

Make sure you install package "Gadfly". If not, you can use "Pkg.add("Gadfly")" to install it.

## API Introduction

```
getFilteredAllSuboptimalAlignment(
    dataPath;
    useAffine = true,
    indelPenalty = -10,
    affineStart = -12,
    affineExtend = -3,
    threshold = 60,
    width = 3,
    plotAllAlignments = false,
    plotGroupAlignments = false
    )
```

### Parameters:
* dataPath: String, the path of two sequences file, data file must be in .fa format
* useAffine: Boolean, set to true if you wish to use Affine penalty strategy. The default value is true
* indelPenalty: Integer, indel penalty for indels, enabled if useAffine = false, the default value is -10
* affineStart: Integer, affine start penalty, enabled if useAffine = true, the defalut value is -12
* affineExtend: Integer, affine extend penalty for every indels after the fist position, enabled if useAffine = true, the default value is -3
* threshold: Integer, threshold for qualified suboptimal alignments
* width: Integer, use to eliminate redundant local alignments, the default value is 3
* plotAllAlignments: Boolean, set to true if you wish to see all local alignments in one graph, the default value is false.
* plotGroupAlignments: Boolean, set to true if you wish to view only repeated suboptimal alignments, the default value is false.

## Usage:

* You can install this package and call getFilteredAllSuboptimalAlignment()
