#!/usr/bin/env python3

import os,sys
import argparse
import time
import scipy.stats as stats
import collections

#######################################
############## ARGUMENTS ##############
#######################################

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument ('--associated','-a',type=str,help='Associated file', required='True')
parser.add_argument ('--window','-w',help='Window',default=5000)
parser.add_argument ('--outfile','-o',type=str, default='./outfile',help ='Output file')

args = parser.parse_args()

##########
## MAIN ##
##########

fileAssociated = open(args.associated,'r')

dictAs = {}
n = 0
for line in fileAssociated:
	n = n + 1
	line = line.strip().split("\t")
	chrom = line[0]
	startcpg = line[1]
	idsnp = line[2]
	startsnp = line[3]
	snp = idsnp+"_"+startsnp
	try:
		asociados = dictAs[chrom,startcpg]
	except:
		asociados = []
	asociados.append(snp)
	dictAs[chrom,startcpg]=asociados

window = int(args.window)
output = {}
for cpg in dictAs:
	bloques = []
	inicio = 0
	end = 0
	chrom = cpg[0]
	startcpg = cpg[1]
	asociados = dictAs[cpg]
	for snp in asociados:
		snp = snp.split("_")
		id = snp[0]
		coordenada = int(snp[1])
		if inicio == 0:
			inicio = coordenada
			fin = coordenada
			ids = id
			continue
		elif (fin+window)>coordenada:
			fin = coordenada
			ids = ids+","+id
		else:
			if inicio!=fin:
				bloque = (inicio,fin,ids)
				bloques.append(bloque)
			inicio = 0
			fin = 0
			ids = ""
	if bloques:
		output[cpg]=bloques

outfile = open(args.outfile,'a')

cabecera = "chrom\tchromStartCpG\tinicioBloque\tfinBloque\tIDS\n"
outfile.write(cabecera)
for cpg in output:
	chrom = str(cpg[0])
	startcpg = str(cpg[1])
	bloques = output[cpg]
	for bloque in bloques:
		inicio = str(bloque[0])
		fin = str(bloque[1])
		ids = str(bloque[2])
		escribir = chrom+"\t"+startcpg+"\t"+inicio+"\t"+fin+"\t"+ids+"\n"
		outfile.write(escribir)
