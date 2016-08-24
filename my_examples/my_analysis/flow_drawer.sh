#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2012-2013 Clemson University.
# 
#  This file was originally written by Michael J. Bojazi and Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file flow_graph.sh
#//! \brief A short shell script to generate abundance graphs.
#//!
#///////////////////////////////////////////////////////////////////////////////

#/bin/bash
#
Nproc=3    # 可同时运行的最大作业数

if [ ! $# -eq 1 ]
  then
    echo -e "\nUsage: $0 basename\n"
    exit
fi

base=$1



#///////////////////////////////////////////////////////////////////////////////
# Edit these parameters.
#///////////////////////////////////////////////////////////////////////////////

INPUT_DIR=dot_flow
OUTPUT_DIR=dot_flow_png

if [ -d "$INPUT_DIR" ]; then
  rm -fr $INPUT_DIR
fi
if [ -d "$OUTPUT_DIR" ]; then
  rm -fr $OUTPUT_DIR
fi

mkdir $INPUT_DIR 
mkdir $OUTPUT_DIR 

function CMD {        # 测试命令, 随机等待几秒钟
    neato -n -Tpng "${inputFilesArray[$1]}" > $OUTPUT_DIR/$(basename "${inputFilesArray[$1]}").png
}

./zone_flow_graph ../network/${base}_output.xml "[z <= 20 and a - z <= 20]" "[product[contains(.,'gamma')] or product[contains(.,'positron')]]" "" "[z > 3 and a - z > 3 and z <= 20 and a - z <= 20]" "all" "linear" dot_flow/flow

#///////////////////////////////////////////////////////////////////////////////
# End of edit.
#///////////////////////////////////////////////////////////////////////////////

ORIGINAL_DIR=`pwd`

mkdir $OUTPUT_DIR

#///////////////////////////////////////////////////////////////////////////////
# Find all the input dot files and put them into an array.
#///////////////////////////////////////////////////////////////////////////////

inputFilesArray=($(find $INPUT_DIR -type f))

#///////////////////////////////////////////////////////////////////////////////
# Now, find the length of the array.
#///////////////////////////////////////////////////////////////////////////////

numfiles=${#inputFilesArray[@]}

#///////////////////////////////////////////////////////////////////////////////
# Loop over all the files in the array.
#///////////////////////////////////////////////////////////////////////////////
PID=() # 记录PID到数组, 检查PID是否存在以确定是否运行完毕
for(( i=0; i<numfiles; )); 
do  
    for((Ijob=0; Ijob<Nproc; Ijob++)); do
        if [[ $i -ge $numfiles ]]; then
            break;
        fi
        if [[ ! "${PID[Ijob]}" ]] || ! kill -0 ${PID[Ijob]} 2> /dev/null; then
            CMD $i $Ijob &
            PID[Ijob]=$!
            i=$((i+1))
        fi
    done
    sleep 0.01
done  
wait

outputFilesArray=("$(ls $OUTPUT_DIR | sort -n -t _ -k 5)")

cd $OUTPUT_DIR

# fmpeg -nostdin -v -1 -i '$(basename "${inputFilesArray[$i]}").png' 'dot.mp4'

#///////////////////////////////////////////////////////////////////////////////
# Use ghostscript to combine all the separate pdf files into a single pdf file.
#///////////////////////////////////////////////////////////////////////////////

# gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${ORIGINAL_DIR}/$5 ${outputFilesArray[@]}
ffmpeg -nostdin -v -1 -i 'flow_%5d.png' "../${base}_FLOW.mp4"

cd "$ORIGINAL_DIR"

#///////////////////////////////////////////////////////////////////////////////
# Remove these lines if you want to retain the individual dot and tex files.
#///////////////////////////////////////////////////////////////////////////////

# rm -rf $INPUT_DIR 
# rm -fr $OUTPUT_DIR
