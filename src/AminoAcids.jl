type AminoAcid
    name
    length
    content
end

type Alignment
    seqX # string
    seqY # string
    startPos   # (x, y)
    endPos
    score
end

type Group
    startCoord
    endCoord
    members
    memberNum
end

aminoAcids = Set()
push!(aminoAcids,
"A",
"R",
"N",
"D",
"B",
"C",
"E",
"Q",
"Z",
"G",
"H",
"I",
"L",
"K",
"M",
"F",
"P",
"S",
"T",
"W",
"Y",
"V")

function loadData(dataPath)
    data = open(dataPath)
    lines = readlines(data)

    aminoAcidlist = []

    for i = 1:2:length(lines)
        lineEles = split(lines[i], ",")
        contend = split(lines[1+i], "")

        for j = contend
            if !checkElement(j, aminoAcids)
                print("Wrong element: ")
                print(j)
                print("\n")
                return
            end
        end

        aminoAcid = AminoAcid(
            split(lineEles[1])[2],
            split(lineEles[2])[2],
            split(lines[1+i], "")
        )
        push!(aminoAcidlist, aminoAcid)
    end

    return aminoAcidlist
end

function checkElement(target, ele)
    for i = ele
        if i == target
            return true
        end
    end
    return false
end



#=
A   ALA
R   ARG
N   ASN
D   ASP
C   CYS
Q   GLN
E   GLU
G   GLY
H   HIS
I   ILE
L   LEU
K   LYS
M   MET
F   PHE
P   PRO
S   SER
T   THR
W   TRP
Y   TYR
V   VAL

A   1
R   2
N   3
D   4
C   5
Q   6
E   7
G   8
H   9
I   10
L   11
K   12
M   13
F   14
P   15
S   16
T   17
W   18
Y   19
V   20

4   -1  -2  -2  0   -1  -1  0   -2  -1  -1  -1  -1  -2  -1  1   0   -3  -2  0
-1  5   0   -2  -3  1   0   -2  0   -3  -2  2   -1  -3  -2  -1  -1  -3  -2  -3
-2  0   6   1   -3  0   0   0   1   -3  -3  0   -2  -3  -2  1   0   -4  -2  -3
-2  -2  1   6   -3  0   2   -1  -1  -3  -4  -1  -3  -3  -1  0   -1  -4  -3  -3
0   -3  -3  -3  9   -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1
-1  1   0   0   -3  5   2   -2  0   -3  -2  1   0   -3  -1  0   -1  -2  -1  -2
-1  0   0   2   -4  2   5   -2  0   -3  -3  1   -2  -3  -1  0   -1  -3  -2  -2
0   -2  0   -1  -3  -2  -2  6   -2  -4  -4  -2  -3  -3  -2  0   -2  -2  -3  -3
-2  0   1   -1  -3  0   0   -2  8   -3  -3  -1  -2  -1  -2  -1  -2  -2  2   -3
-1  -3  -3  -3  -1  -3  -3  -4  -3  4   2   -3  1   0   -3  -2  -1  -3  -1  3
-1  -2  -3  -4  -1  -2  -3  -4  -3  2   4   -2  2   0   -3  -2  -1  -2  -1  1
-1  2   0   -1  -3  1   1   -2  -1  -3  -2  5   -1  -3  -1  0   -1  -3  -2  -2
-1  -1  -2  -3  -1  0   -2  -3  -2  1   2   -1  5   0   -2  -1  -1  -1  -1  -1
-2  -3  -3  -3  -2  -3  -3  -3  -1  0   0   -3  0   6   -4  -2  -2  1   3   1
-1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4  7   -1  -1  -4  -3  -2
1   -1  1   0   -1  0   0   0   -1  -2  -2  0   -1  -2  -1  4   1   -3  -2  -2
0   -1  0   -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1  1   5   -2  -2  0
-3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1  1   -4  -3  -2  11  2   -3
-2  -2  -2  -3  -2  -1  -2  -3  2   -1  -1  -2  -1  3   -3  -2  -2  2   7   -1
0   -3  -3  -3  -1  -2  -2  -3  -3  3   1   -2  -1  1   -2  -2  0   -3  -1  4

=#

using DataFrames

function getBlosum()
    dataPath = joinpath(@__DIR__, "..", "data")
    blosumPath = joinpath(dataPath, "ezBlosum.csv")
    blosum = readtable(blosumPath, separator = ',')
    return blosum
end

blosum = getBlosum()

blosumIndex = Dict(
    "A" => 1,
    "R" => 2,
    "N" => 3,
    "D" => 4,
    "C" => 5,
    "Q" => 6,
    "E" => 7,
    "G" => 8,
    "H" => 9,
    "I" => 10,
    "L" => 11,
    "K" => 12,
    "M" => 13,
    "F" => 14,
    "P" => 15,
    "S" => 16,
    "T" => 17,
    "W" => 18,
    "Y" => 19,
    "V" => 20
)

function getMappingScore(col, row, blosum)
    col = uppercase(col)
    row = uppercase(row)

    col = blosumIndex[col]
    row = blosumIndex[row]

    return blosum[col][row]
end

#=
This function searches for the maximum score left in the score matrix. If
it is >= the score threshold, then it will report the alignment corresponding
to it and return true, else false is returned.
@param matrix Score matrix
@param path Path matrix
@param xLen Length of seqX + 1
@param yLen Length of seqY + 1
@param alignmentList List of alignments to which we will append new alignment
@param threshold Lowest allowable score
@return True if alignment has sufficient score and is successfully added to
        list, else false
=#
function getAlignment(matrix, path, xLen, yLen, alignmentList, threshold, width, xAminoAcid, yAminoAcid)
    maxScore = 0
    col = 0
    row = 0
    ## First get the maxScore and its coordinates
    for r = 2:yLen
        for c = 2:xLen
            if matrix[r, c] > maxScore
                maxScore = matrix[r, c]
                col = c
                row = r
            end
        end
    end
    ## Check if maxScore >= threshold
    if maxScore < threshold
        return false ## score too low
    end
    ## Now retrieve the path
    seqX = ""
    seqY = ""
    x1 = 0
    x2 = col
    lastX = col
    y1 = 0
    y2 = row
    lastY = row

    #zero out the region around the maxScore
    for w in 0:width
        if row - w >= 2 && col - w >= 2
            matrix[row - w, col - w] = 0
        end
        if row - w >= 2 && col + w <= xLen
            matrix[row - w, col + w] = 0
        end
        if row + w <= yLen && col - w >= 2
            matrix[row + w, col - w] = 0
        end
        if row + w <= yLen && col + w <= xLen
            matrix[row + w, col + w] = 0
        end
    end
    pathVal = path[row, col] ## init path location
    while pathVal != 0
        ## First zero out the score element and nearby left/right region
        matrix[row, col] = 0
        for w in 1:width
            if col - w >= 2
                matrix[row, col - w] = 0
            end
            if col + w <= xLen
                matrix[row, col + w] = 0
            end
        end
        lastX = row
        lastY = col
        if pathVal == 1  # match/mismatch
            row -= 1
            col -= 1
            seqX = string(xAminoAcid.content[col], seqX)
            seqY = string(yAminoAcid.content[row], seqY)
        elseif pathVal == 2
            row -= 1
            seqX = string("-", seqX)
            seqY = string(yAminoAcid.content[row], seqY)
        else ## pathVal == 3
            col -= 1
            seqX = string(xAminoAcid.content[col], seqX)
            seqY = string("-", seqY)
        end
        pathVal = path[row, col]
    end
    x1 = lastX
    y1 = lastY
    ## Add this alignment to the list of alignments
    alignment = Alignment(
        seqX,
        seqY,
        (x1, y1),
        (y2, x2),
        maxScore
    )
    push!(alignmentList, alignment)
    return true
end

#=
This function returns a filtered list of alignments. If alignments are
overlapping (i.e. have the same start position), then only the one with
the highest score will be added.
@param alignmentList List of alignments
@return filtered list of non-overlapping alignments
=#
function filterAlignments(alignmentList)
    # sort the alignments first
    sort!(alignmentList, by = x -> x.startPos)
    filtList = []
    numFilt = 1
    i = length(alignmentList) - 1
    # push / pop the very first alignment
    push!(filtList, pop!(alignmentList))
    while i > 0
        if alignmentList[i].startPos == filtList[numFilt].startPos
            if alignmentList[i].score > filtList[numFilt].score
                pop!(filtList)
                push!(filtList,pop!(alignmentList))
            else
                pop!(alignmentList)
            end
        else # new start position
            push!(filtList,pop!(alignmentList))
            numFilt += 1
        end
        i -= 1
    end
    return filtList
end

#=
 = alignmentList is a list of Alignment object
 =#

function mergeAlignment(alignmentList)
    seqXGroups = [] # list of group x object
    seqYGroups = [] # list of group y object
    for alignment = alignmentList
        pushX = false
        for xGroup in seqXGroups
            if isCoordSimilar((xGroup.startCoord, xGroup.endCoord), (alignment.startPos[2], alignment.endPos[2]))
                push!(xGroup.members, alignment)
                xGroup.startCoord = (xGroup.startCoord * xGroup.memberNum + alignment.startPos[2]) / (xGroup.memberNum + 1)
                xGroup.endCoord = (xGroup.endCoord * xGroup.memberNum + alignment.endPos[2]) / (xGroup.memberNum + 1)
                xGroup.memberNum += 1
                pushX = true
            end
        end

        if pushX == false
            newSet = Set()
            push!(newSet, alignment)
            push!(
                seqXGroups,
                Group(
                    alignment.startPos[2],
                    alignment.endPos[2],
                    newSet,
                    1
                )
            )
        end

        pushY = false
        for yGroup in seqYGroups
            if isCoordSimilar((yGroup.startCoord, yGroup.endCoord), (alignment.startPos[1], alignment.endPos[1]))
                push!(yGroup.members, alignment)
                yGroup.startCoord = (yGroup.startCoord * yGroup.memberNum + alignment.startPos[1]) / (yGroup.memberNum + 1)
                yGroup.endCoord = (yGroup.endCoord * yGroup.memberNum + alignment.endPos[1]) / (yGroup.memberNum + 1)
                yGroup.memberNum += 1
                pushY = true
            end
        end

        if pushY == false
            newSet = Set()
            push!(newSet, alignment)
            push!(
                seqYGroups,
                Group(
                    alignment.startPos[1],
                    alignment.endPos[1],
                    newSet,
                    1
                )
            )
        end
    end
    return (seqXGroups, seqYGroups)
end

#=
    xCoord1 = (x1, x2) / (y1, y2)
    xCoord2 = (x1, x2) / (y1, y2)

    e.g.
    |xCoord1(x1) - xCoord2(x1)| <= distance (usually be 1)
    |xCoord1(x2) - xCoord2(x2)| <= distance (usually be 1)
=#
distance = 2
function isCoordSimilar(xCoord1, xCoord2)
    if (abs(xCoord1[1] - xCoord2[1]) <= distance) && (abs(xCoord1[2] - xCoord2[2]) <= distance)
        return true
    else
        return false
    end
end

function plotAlignmentProfile(alignmentListFilt, path, scoreMatrix)
    mark = 0
    for a in alignmentListFilt
        mark += 1
        row = a.endPos[1]
        col = a.endPos[2]
        score = []
        while path[row, col] != 0
            unshift!(score, scoreMatrix[row, col])
            if path[row, col] == 1
                row -= 1
                col -= 1
            elseif path[row, col] == 2 # indel from left
                col -= 1 
            elseif path[row, col] == 3 # indel from right
                row -= 1
            end
        end
        graphTitle = string("Alignment", dec(mark), " Score Profile")
        
        alignmentScore = plot(x = collect(1:length(score)), y = score, Guide.xlabel("Alignment Position"), Guide.ylabel("Score"), Guide.title(graphTitle))
        fileName = string("Alignment", dec(mark))
        fileName = string(fileName, "ScoreProfile")
        fileName = string(fileName, ".svg")
        draw(SVG(fileName, 3inch, 3inch), alignmentScore)
    end
end





