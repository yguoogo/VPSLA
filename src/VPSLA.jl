module VPSLA

# package code goes here

include("AminoAcids.jl")
using Gadfly
#### Parameters ####
#=indelPenalty = -3
affineStart = -100
affineExtend = -1
threshold = 30
width = 3 =#
####################

export getFilteredAllSuboptimalAlignment

function plotAllSuboptimalAlignment(alignmentListFilt, xAminoAcid, yAminoAcid, path)
    xAxis = copy(xAminoAcid.content)
    unshift!(xAxis, "0")
    yAxis = copy(yAminoAcid.content)
    unshift!(yAxis, "0")

    xLen = length(xAminoAcid.content)+1
    yLen = length(yAminoAcid.content)+1

    allShowPathMatrix = getAllShowPathMatrix(alignmentListFilt, xLen, yLen, path)
    draw(SVG("allPath.svg", 20inch, 20inch),
    spy(allShowPathMatrix, Scale.y_discrete(labels = i->yAxis[i]), Scale.x_discrete(labels = i->xAxis[i])))
end

function plotSuboptimalAlignmentGroup(res, xAminoAcid, yAminoAcid, path)
    xAxis = copy(xAminoAcid.content)
    unshift!(xAxis, "0")
    yAxis = copy(yAminoAcid.content)
    unshift!(yAxis, "0")

    xLen = length(xAminoAcid.content)+1
    yLen = length(yAminoAcid.content)+1

    for group = res[1]
        pathMark = 0
        if group.memberNum > 1
            pathMark += 1
            showPathMatrix = getGroupShowPathMatrix(group.members, pathMark, xLen, yLen, path)
            fileName = string("GroupX", dec(pathMark))
            fileName = string(fileName, ".svg")
            draw(SVG(fileName, 20inch, 20inch),
            spy(showPathMatrix, Scale.y_discrete(labels = i->yAxis[i]), Scale.x_discrete(labels = i->xAxis[i])))
        end
    end

    for group = res[2]
        pathMark = 0
        if group.memberNum > 1
            pathMark += 1
            showPathMatrix = getGroupShowPathMatrix(group.members, pathMark, xLen, yLen, path)
            fileName = string("GroupY", dec(pathMark))
            fileName = string(fileName, ".svg")
            draw(SVG(fileName, 20inch, 20inch),
            spy(showPathMatrix, Scale.y_discrete(labels = i->yAxis[i]), Scale.x_discrete(labels = i->xAxis[i])))
        end
    end

end

function getGroupShowPathMatrix(alignmentList, pathMark, xLen, yLen, path)

    showPath = zeros(yLen, xLen)

    for a in alignmentList
        row = a.endPos[1]
        col = a.endPos[2]

        while path[row, col] != 0
            showPath[row, col] = pathMark

            if path[row, col] == 1
                row -= 1
                col -= 1
            elseif path[row, col] == 2 # indel from left
                col -= 1
            elseif path[row, col] == 3 # indel from right
                row -= 1
            end
        end
    end

    return showPath
 end

 function getAllShowPathMatrix(alignmentList, xLen, yLen, path)
    pathMark = 0
    showPath = zeros(yLen, xLen)

    for a in alignmentList
        pathMark += 1
        row = a.endPos[1]
        col = a.endPos[2]

        while path[row, col] != 0
            showPath[row, col] = pathMark

            if path[row, col] == 1
                row -= 1
                col -= 1
            elseif path[row, col] == 2 # indel from left
                col -= 1
            elseif path[row, col] == 3 # indel from up
                row -= 1
            end
        end
    end
    return showPath
 end

 function reportGroup(res)
    for group = res[1]
        pathMark = 0
        if group.memberNum > 1
            pathMark += 1
            println("X axis Alignment Group$(pathMark)")
            for member in group.members
                println("Score = $(member.score), start = $(member.startPos), end = $(member.endPos):")
                println(member.seqX)
                println(member.seqY)
                println("start position = $(member.startPos)")
                println("end position = $(member.endPos)")
                println()
            end
        end
    end

    for group = res[2]
        pathMark = 0
        if group.memberNum > 1
            pathMark += 1
            println("Y axis Alignment Group$(pathMark):")
            for member in group.members
                println("Score = $(member.score), start = $(member.startPos), end = $(member.endPos):")
                println(member.seqX)
                println(member.seqY)
                println("start position = $(member.startPos)")
                println("end position = $(member.endPos)")
                println()
            end
        end
    end
end

function printAllAlignments( alignmentListFilt )
    println("---------------------------------------------------------")
    println("Printing all alignments...")
    println()
    for a in alignmentListFilt
        println("Score = $(a.score), start = $(a.startPos), end = $(a.endPos):")
        println(a.seqX)
        println(a.seqY)
        println("start position = $(a.startPos)")
        println("end position = $(a.endPos)")
        println()
    end
    println("---------------------------------------------------------")
end

function getFilteredAllSuboptimalAlignment(dataPath; indelPenalty = -10, affineStart = -12, affineExtend = -3, threshold = 30, width = 3, useAffine = true, plotAllAlignments = false, plotGroupAlignments = false)

    list = loadData(dataPath)
    blosum = getBlosum()
    xAminoAcid = list[1]
    yAminoAcid = list[2]
    xLen = length(xAminoAcid.content)+1
    yLen = length(yAminoAcid.content)+1
    matrix = zeros(yLen, xLen)

    path = zeros(yLen, xLen)

    for row = 2:yLen
        for col = 2:xLen

            #=
             = In path matrix, 1 stands for match/mismatch,
             = 2 stands for indel from left
             = 3 stands for indel from up
             =#

            # if we chose the path with match/mismatch
            a = matrix[row-1, col-1] + getMappingScore(yAminoAcid.content[row-1], xAminoAcid.content[col-1], blosum)
            # calculate indel

            b = 0
            if useAffine
                if path[row-1, col] == 2 # There already indels from left at previous left cell
                    b = matrix[row-1, col] + affineExtend # 2: indel from left
                else
                    b = matrix[row-1, col] + affineStart
                end
            else
                b = matrix[row-1, col] + indelPenalty
            end

            c = 0
            if useAffine
                if path[row, col-1] == 3 # There already indels from up at previous up cell
                    c = matrix[row, col-1] + affineExtend
                else
                    c = matrix[row, col-1] + affineStart # 3: indel from up
                end
            else
                c = matrix[row, col-1] + indelPenalty
            end

            max = 0
            if (a > max)
                max = a
                path[row, col] = 1
            end
            if (b > max)
                max = b
                path[row, col] = 2
            end
            if (c > max)
                max = c
                path[row, col] = 3
            end
            if (0 > max)
                max = 0
                path[row, col] = 0
            end
            matrix[row, col] = max
        end
    end

    scoreMatrix = copy(matrix)
    ## Get all of the alignments above the threshold score in a list.
    alignmentList = []
    while getAlignment(matrix, path, xLen, yLen, alignmentList, threshold, width, xAminoAcid, yAminoAcid)
    end

    if length(alignmentList) == 0
        println("No local alignments. You can try other parameters like lower the threshold")
        return
    end

    ## Filter out overlapping suboptimal alignments
    alignmentListFilt = filterAlignments(alignmentList)
    res = mergeAlignment(alignmentListFilt)
    reportGroup(res)

    if plotAllAlignments
        plotAllSuboptimalAlignment(alignmentListFilt, xAminoAcid, yAminoAcid, path)
    end

    if plotGroupAlignments
        res = mergeAlignment(alignmentListFilt)
        reportGroup(res)
        plotSuboptimalAlignmentGroup(res, xAminoAcid, yAminoAcid, path)
    end
    plotAlignmentProfile(alignmentListFilt, path, scoreMatrix)
    ## Print all of the alignments
    printAllAlignments( alignmentListFilt ) 
    return alignmentListFilt
end

end # module
