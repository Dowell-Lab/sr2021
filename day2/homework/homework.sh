#!/usr/bin/env bash

# This homework is going to walk you through a few different exercises
# using some of the command line tools we talked about in class.

# Please don't be afraid to open up the files and look at them!
# It's a great way to figure out what's actually going on in the files.

################################################################################
# Question 1:
#
# We're going to start by using cat and grep to manipulate some files
# containing several lines. Take a look at the files
# heart_shaped_box.txt and star_shaped_box.txt in the folder boxes.
################################################################################
# First, let's combine the heart shaped and star shaped box files into
# a single file called combined_box.txt using the cat command:
# YOUR CODE HERE
cat heart_shaped_box.txt star_shaped_box.txt > combined_box.txt

# Next, we'll Write a new file called polka_dot_number.txt that
# contains the number of polka dot balls from the combined box.

# One way we can do this is to use grep to find lines that have the
# word polka_dot in them, and then send those to a file. Use grep to
# find the lines that have polka_dot on them, and output the result to
# a file called polka_dot.txt in the boxes directory.
# YOUR CODE HERE
grep "polka_dot" combined_box.txt > polka_dot.txt
# When we're working with real data, we don't necessarily want to have
# to manually count lines every time we search for something. Instead,
# use the wc (word count) utility to count the number of lines our
# search returns using wc -l (-l = count lines). Output the result to
# a file called polka_dot_number.txt in the boxes directory.
# YOUR CODE HERE
wc -l polka_dot.txt > polka_dot_number.txt

################################################################################
# Question 2
#
# In the next part of this homework, we are going to look at some
# files in the same format as the bedgraph files that are actually
# used in the analysis of sequencing data.
################################################################################

# First, we want to make an output directory to hold the files we
# produce. For this question we will be working in the 'genes'
# directory, so let's make an output directory called out inside of
# the genes directory using mkdir:
# YOUR CODE HERE
mkdir -p genes/out

# Note:
#
# For this question, We are looking at tab separated text files with
# the format:
#chromosome start stop gene_name exon #

# Let's start by finding all of the exons for chromosome 1 in humans
# and outputting them in a file called human_chr1_exons.txt inside the
# output directory we just made:
# YOUR CODE HERE
grep "chr1" genes/

# Next, let's Write out a file that contains the total number of exons
# in chromosome 1 for mice, putting it in a file called
# mouse_chr1_exon_count.txt inside the output directory we just made:
# YOUR CODE HERE

# Are there the same number of genes in chromosome 1 for mice and
# humans? Using the human_genes.txt and mouse_genes.txt files,
# see if this is the case.
# HINT (Which exon is always present?)
# YOUR CODE HERE
# YOUR CODE HERE

# You notice there is alternative splicing in the sonic_hedgehog gene
# in mouse compared to humans. Write a file containing only the stop
# coordinate of the last exon for this gene in mice. Do the same for
# human. Are these coordinates the same?
# HINT (Look at the last row.)
# HINT (It might be useful to use the awk utility to select a specific
# column)
# YOUR CODE HERE
# YOUR CODE HERE

# You notice that worm researchers have put unnecessary comments
# starting with "#" into their genes file. Use grep to take out the
# comments. Put the output in a file called worm_genes_no_comments.txt
# in the output directory we made previously.
# YOUR CODE HERE

# You notice that for some reason the gene files you have use a
# different chromosome system than everyone else. Add "chr" in front
# of the numbers to the first column of the file you just removed
# comments from so that it matches your other files. Put the result in
# a file called worm_genes_fixed.txt in the output directory we already
# made.
# YOUR CODE HERE

# You also noticed that the files you have are highly unorganized and
# aren't sorted by chromosome or start coordinate. Fix this by sorting
# your files using the sort command.
# HINT (sort -k1,1 -kn2,2)
# YOUR CODE HERE

################################################################################
# Advanced Questions:
#
# These are a few more advanced questions showing off a few different
# more advanced features of the command line.
################################################################################

# Are there any genes in chromosome 2 that have the same name in mice
# and worms?

# There are many ways you could do this. One way would be to get files
# with only chr2 gene names, then use a while loop to iterate
# through one file and grep the gene name from the other. If the gene
# name is found, then append it to the new file.

# To start, select only genes in chromosome 2:
grep chr2 genes/mouse_genes.txt | grep exon1 | awk '{ print $4 }' > genes/out/mouse_genes_chr2.txt
grep chr2 genes/out/worm_genes_final.txt | grep exon1 | awk '{ print $4 }' > genes/out/worm_genes_chr2.txt

# The first way we can solve the problem above is to use a while loop
# as described above.

# Create an empty file so you don't keep appending if you run the bash
# loop again
echo '' > genes/out/mouse_worm_genes_chr2.txt
# Make variables for the files we're going to use in the loop so we
# don't have to keep writing the full filename.
file1=genes/out/mouse_genes_chr2.txt
file2=genes/out/worm_genes_chr2.txt
# Run the loop. p is the current line we read in from file1. In the
# loop, we search to see if that same gene name is in the worm genes
# file (file2).
while read -r p; do
		echo grepping gene: "$p"
		grep "$p" "$file2" >> genes/out/mouse_worm_genes_chr2.txt
done < "$file1" # < "$file1" feeds file1 into the read function in our
								# loop line by line

# If we don't want to use a loop, we can also solve this problem using
# just pipes. We do this using the following steps:
# 1. Combine both gene lists into a single file.
# 2. Sort the combined gene list so that duplicate genes are next to
#    each other
# 3. Count the unique elements in the sorted, combined list. If a gene
#    is duplicated, it will be counted more than once.
# 4. Using awk, check if the first column (the count) is greater than 1.
#    If it is, the gene must have been in both files so we print out the
#    gene name (the second column).
# 5. Write the unique gene names out to a file.
cat genes/out/mouse_genes_chr2.txt genes/out/worm_genes_chr2.txt | sort | \
		uniq -c | awk '{if ($1 > 1) print $2}' \
									> genes/out/mouse_worm_genes_chr2_zach.txt
# Although this method is a little bit more complex than the method
# using a while loop, it is also much faster (especially when your
# gene files get large). In the while loop example, we look through
# the entire worm genes file for every single mouse gene we have. If
# both files have 1000 genes in them, this means we have to look at
# 1000*1000=1000000 (a million!) lines. In the second method using
# pipes, we only have to read in the files once and look over the
# entire list once for each step in the pipeline.
