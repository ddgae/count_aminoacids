#!/usr/bin/env bash
#
# Author: David Gae @ UCSF
# Citation: Boyken et al. 2016 De Novo Design of protein homo-oligomers with modular hydrogen-bond network-mediated specificity
#           Wu et al. 2000 Coiled Trigger Motifs in the 1B and 2B Rod Domain Segments Are Required for the Stability of Keratin Intermediate Filaments 
#	    Lumb et al. 1995 A Buried Polar Interaction Imparts Structural Uniqueness in a Designed Heterodimeric Coiled Coilf
#           Kuster et al. 2015 High-Resolution Crystal Structures of Protein Helices Reconciled with Three-Centered Hydrogen Bonds and Multipole Electrostatics
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# HYDROGEN PATTERN ANALYSIS, and user specific acceptor group and, its interaction to all possible donors in the
# protein and a sample ligand.
# determine Hydrogen Bond pattern
# usage: python hbond_pattern.py "input.pdb" "acceptor"
# LIMITATION: only works for one or two chain protein such as an abeta protein and one substrate. still need improvement
# on creating option commands if more chains are present or atomtype is different then suggested.
import sys,os

if sys.argv[0]: print("Example: python hbond_pattern.py template.pdb hbond_acceptor A B")

#if sys.argv[0]:
#   sys.exit()

import math
import sys
import numpy as np
import re

#ACCEPTOR-DONOR ASN (4)
#ASNOD1_ARGNE
#ASNOD1_ARGNH1
#ASNOD1_ARGNH2
#ASNOD1_ASNND2
#ASNOD1_GLNNE2
#ASNOD1_HISND1
#ASNOD1_HISNE2
#ASNOD1_LYSNZ
#ASNOD1_SEROG
#ASNOD1_THROG1
#ASNOD1_TRPNE1
#ASNOD1_TYROH

#ACCEPTOR-DONOR GLN
#GLNOE1_ARGNE
#GLNOE1_ARGNH1
#GLNOE1_ARGNH2
#GLNOE1_ASNND2
#GLNOE1_GLNNE2
#GLNOE1_HISND1
#GLNOE1_HISNE2
#GLNOE1_LYSNZ
#GLNOE1_SEROG
#GLNOE1_THROG1
#GLNOE1_TRPNE1
#GLNOE1_TYROH

#DONE

#ACCEPTOR-DONOR GLU
#GLUOE1_ARGNE
#GLUOE2_ARGNE
#GLUOE1_ARGNH1
#GLUOE2_ARGNH1
#GLUOE1_ARGNH2
#GLUOE2_ARGNH2
#GLUOE1_ASNND2
#GLUOE2_ASNND2
#GLUOE1_GLNNE2
#GLUOE2_GLNNE2
#GLUOE1_HISND1
#GLUOE1_HISNE2
#GLUOE2_HISND1
#GLUOE2_HISNE2
#GLUOE1_LYSNZ
#GLUOE2_LYSNZ
#GLUOE1_SEROG
#GLUOE2_SEROG
#GLUOE1_THROG1
#GLUOE2_THROG1
#GLUOE1_TRPNE1
#GLUOE2_TRPNE1
#GLUOE1_TYROH
#GLUOE2_TYROH

#ACCEPTOR-DONOR ASP
#ASPOD1_ARGNE
#ASPOD2_ARGNE
#ASPOD1_ARGNH1
#ASPOD2_ARGNH1
#ASPOD1_ARGNH2
#ASPOD2_ARGNH2
#ASPOD1_ASNND2
#ASPOD2_ASNND2
#ASPOD1_GLNNE2
#ASPOD2_GLNNE2
#ASPOD1_HISND1
#ASPOD1_HISNE2
#ASPOD2_HISND1
#ASPOD2_HISNE2
#ASPOD1_LYSNZ
#ASPOD2_LYSNZ
#ASPOD1_SEROG
#ASPOD2_SEROG
#ASPOD1_THROG1
#ASPOD2_THROG1
#ASPOD1_TRPNE1
#ASPOD2_TRPNE1
#ASPOD1_TYROH
#ASPOD2_TYROH

#ACCEPTOR-DONOR HSD
#HISND1_ARGNE
#HISNE2_ARGNE
#HISND1_ARGNH1
#HISNE2_ARGNH1
#HISND1_ARGNH2
#HISNE2_ARGNH2
#HISND1_ASNND2
#HISNE2_ASNND2
#HISND1_GLNNE2
#HISNE2_GLNNE2
#HISND1_HISND1
#HISND1_HISNE2
#HISNE2_HISND1
#HISNE2_HISNE2
#HISND1_LYSNZ
#HISNE2_LYSNZ
#HISND1_SEROG
#HISNE2_SEROG
#HISND1_THROG1
#HISNE2_THROG1
#HISND1_TRPNE1
#HISNE2_TRPNE1
#HISND1_TYROH
#HISNE2_TYROH

#ACCEPTOR-DONOR SER
#SEROG_ARGNE
#SEROG_ARGNH1
#SEROG_ARGNH2
#SEROG_ASNND2
#SEROG_GLNNE2
#SEROG_HISND1
#SEROG_HISNE2
#SEROG_LYSNZ
#SEROG_SEROG
#SEROG_THROG1
#SEROG_TRPNE1
#SEROG_TYROH

#ACCEPTOR-DONOR THR
#THROG1_ARGNE
#THROG1_ARGNH1
#THROG1_ARGNH2
#THROG1_ASNND2
#THROG1_GLNNE2
#THROG1_HISND1
#THROG1_HISNE2
#THROG1_LYSNZ
#THROG1_SEROG
#THROG1_THROG1
#THROG1_TRPNE1
#THROG1_TYROH

#ACCEPTOR_DONOR TYR
#TYROH_ARGNE
#TYROH_ARGNH1
#TYROH_ARGNH2
#TYROH_ASNND2
#TYROH_GLNNE2
#TYROH_HISND1
#TYROH_HISNE2
#TYROH_LYSNZ
#TYROH_SEROH
#TYROH_THROG1
#TYROH_TRPNE1
#TYROH_TYROH

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

Xvalue1 = []
Yvalue1 = []
Xvalue2 = []
Yvalue2 = []
Xvalue3 = []
Yvalue3 = []
Xvalue4 = []
Yvalue4 = []
Xvalue5 = []
Yvalue5 = []
Xvalue6 = []
Yvalue6 = []
Xvalue7 = []
Yvalue7 = []
Xvalue8 = []
Yvalue8 = []
Xvalue9 = []
Yvalue9 = []
Xvalue10 = []
Yvalue10 = []
Xvalue11 = []
Yvalue11 = []
Xvalue12 = []
Yvalue12 = []
Xvalue13 = []
Yvalue13 = []
Xvalue14 = []
Yvalue14 = []
Xvalue15 = []
Yvalue15 = []
Xvalue16 = []
Yvalue16 = []
Xvalue17 = []
Yvalue17 = []

dist1 = []
dist2 =[]
dist3 = []
dist4 = []
dist5 = []
dist6= []
dist7 = []
dist8 = []
dist9 = []
dist10 = []
dist11 = []
dist12 = []
dist13 = []
dist14 =[]
dist15 = []
dist16 = []
dist17= []
dist18 = []
dist19 = []
dist20 = []
dist21 = []
dist22 = []
dist23 = []
dist24 = []
dist25 = []
dist26 = []
dist27 = []
dist28 = []
dist29 = []
dist30 = []
dist31 = []
dist32 = []
dist33 = []
dist34 = []
dist35 = []
dist36 = []
dist37 = []
dist38 = []
dist39 = []
dist40 = []
dist41 = []
dist42 = []
dist43 = []
dist44 = []
dist45 = []
dist46 = []
dist47 = []
dist48 = []
dist49 = []
dist50 = []
dist51 = []
dist52 = []
dist53 = []
dist54 = []
dist55 = []
dist56 = []
dist57 = []
dist58 = []
dist59 = []
dist60 = []
dist61 = []
dist62 = []
dist63 = []
dist64 = []
dist65 = []
dist66 = []
dist67 = []
dist68 = []
dist69 = []
dist70 = []
dist71 = []
dist72 = []
dist73 = []
dist74 = []
dist75 = []
dist76 = []
dist77 = []
dist78 = []
dist79 = []
dist80 = []
dist81 = []
dist82 = []
dist83 = []
dist84 = []
dist85 = []
dist86 = []
dist87 = []
dist88 = []
dist89 = []
dist90 = []
dist91 = []
dist92 = []
dist93 = []
dist94 = []
dist95 = []
dist96 = []
dist97 = []
dist98 = []
dist99 = []
dist100 = []
dist101 = []
dist102 = []
dist103 = []
dist104 = []
dist105 = []
dist106 = []
dist107 = []
dist108 = []
dist109 = []
dist110 = []
dist111 = []
dist112 = []
dist113 = []
dist114 = []
dist115 = []
dist116 = []
dist117 = []
dist118 = []
dist119 = []
dist120 = []
dist121 = []
dist122 = []
dist123 = []
dist124 = []
dist125 = []
dist126 = []
dist127 = []
dist128 = []
dist129 = []
dist130 = []
dist131 = []
dist132 = []
dist133 = []
dist134 = []
dist135 = []
dist136 = []
dist137 = []
dist138 = []
dist139 = []
dist140 = []
dist141 = []
dist142 = []
dist142a = []
dist143 = []
dist144 = []
dist145 = []
dist146 = []
dist147 = []
dist148 = []
dist149 = []
dist150 = []
dist151 = []
dist152 = []
dist153 = []
dist154 = []
dist155 = []
dist156 = []
dist157 = []
dist158 = []
dist159 = []
dist160 = []
dist161 = []
dist162 = []
dist163 = []
dist164 = []
dist165 = []
dist166 = []
dist167 = []
dist168 = []
dist169 = []
dist170 = []
dist171 = []
dist172 = []
dist173 = []
dist174 = []
dist175 = []
dist176 = []
dist177 = []
dist178 = []
dist179 = []
dist180 = []
dist181 = []
dist182 = []
dist183 = []
dist184 = []
dist185 = []
dist186 = []
dist187 = []
dist188 = []
dist189 = []
dist190 = []
dist191 = []
dist192 = []
dist193 = []
dist194 = []
dist195 = []
dist196 = []
dist197 = []
dist198 = []
dist199 = []
dist200 = []
dist201 = []
dist202 = []
dist203 = []
dist204 = []
dist205 = []
dist206 = []
dist207 = []
dist208 = []
dist209 = []
dist210 = []
dist211 = []
dist212 = []
dist213 = []
dist214 = []
dist215 = []
dist216 = []
dist217 = []
dist218 = []
dist219 = []
dist220 = []
dist221 = []
dist222 = []
dist223 = []
dist224 = []
dist225 = []
dist226 = []
dist227 = []
dist228 = []
dist229 = []
dist230 = []
dist231 = []
dist232 = []
dist233 = []
dist234 = []
dist235 = []
dist236 = []
dist237 = []
dist238 = []
dist239 = []
dist240 = []
dist241 = []
dist242 = []
dist243 = []
dist244 = []
dist245 = []
dist246 = []
dist247 = []
dist248 = []
dist249 = []
dist250 = []
dist251 = []
dist252 = []
dist253 = []
dist254 = []
dist255 = []
dist256 = []
dist257 = []
dist258 = []
dist259 = []
dist260 = []
dist261 = []
dist262 = []
dist263 = []
dist264 = []
dist265 = []
dist266 = []
dist267 = []
dist268 = []
dist269 = []
dist270 = []
dist271 = []
dist272 = []
dist273 = []
dist274 = []
dist275 = []
dist276 = []
dist277 = []
dist278 = []
dist279 = []
dist280 = []
dist281 = []
dist282 = []
dist283 = []
dist284 = []
dist285 = []
dist286 = []
dist287 = []
dist288 = []
dist289 = []
dist290 = []
dist291 = []
dist292 = []
dist293 = []
dist294 = []
dist295 = []
dist296 = []
dist297 = []
dist298 = []
dist299 = []
dist300 = []
dist301 = []
dist302 = []
dist303 = []
dist304 = []
dist305 = []
dist306 = []
dist307 = []
dist308 = []
dist309 = []
dist310 = []
dist311 = []
dist312 = []
dist313 = []
dist314 = []
dist315 = []
dist316 = []
dist317 = []
dist318 = []
dist319 = []
dist320 = []
dist321 = []
dist322 = []
dist323 = []
dist324 = []
dist325 = []
dist326 = []
dist327 = []
dist328 = []
dist329 = []
dist330 = []
dist331 = []
dist332 = []
dist333 = []
dist334 = []
dist335 = []
dist336 = []
dist337 = []
dist338 = []
dist339 = []
dist340 = []
dist341 = []
dist342 = []
dist343 = []
dist344 = []
dist345 = []
dist346 = []
dist347 = []
dist348 = []
dist349 = []
dist350 = []
dist351 = []
dist352 = []
dist353 = []
dist354 = []
dist355 = []
dist356 = []
dist357 = []
dist358 = []
dist359 = []
dist360 = []
dist361 = []
dist362 = []
dist363 = []
dist364 = []
dist365 = []
dist366 = []
dist367 = []
dist368 = []
dist369 = []
dist370 = []
dist371 = []
dist372 = []
dist373 = []
dist374 = []
dist375 = []
dist376 = []
dist377 = []
dist378 = []
dist379 = []
dist380 = []
dist381 = []
dist382 = []
dist383 = []
dist384 = []
dist385 = []
dist386 = []
dist387 = []
dist388 = []
dist389 = []
dist390 = []
dist391 = []
dist392 = []
dist393 = []
dist394 = []
dist395 = []
dist396 = []
dist397 = []
dist398 = []
dist399 = []
dist400 = []
dist401 = []
dist402 = []
dist403 = []
dist404 = []
dist405 = []
dist406 = []
dist407 = []
dist408 = []
dist409 = []
dist410 = []
dist411 = []
dist412 = []
dist413 = []
dist414 = []
dist415 = []
dist416 = []
dist417 = []
dist418 = []
dist419 = []
dist420 = []
dist421 = []
dist422 = []
dist423 = []
dist424 = []
dist425 = []
dist426 = []
dist427 = []
dist428 = []
dist429 = []
dist430 = []
dist431 = []
dist432 = []
dist433 = []
dist434 = []
dist435 = []
dist436 = []
dist437 = []
dist438 = []
dist439 = []
dist440 = []
dist441 = []
dist442 = []
dist443 = []
dist444 = []
dist445 = []
dist446 = []
dist447 = []
dist448 = []
dist449 = []
dist450 = []
dist451 = []
dist452 = []
dist453 = []
dist454 = []
dist455 = []
dist456 = []
dist457 = []
dist458 = []
dist459 = []
dist460 = []
dist461 = []
dist462 = []
dist463 = []
dist464 = []
dist465 = []
dist466 = []
dist467 = []
dist468 = []
dist469 = []
dist470 = []
dist471 = []
dist472 = []
dist473 = []
dist474 = []
dist475 = []
dist476 = []
dist477 = []
dist478 = []
dist479 = []
dist480 = []
dist481 = []
dist482 = []   
dist483 = []
dist484 = []
dist485 = []
dist486 = []
dist487 = []
dist488 = []
dist489 = []
dist490 = []
dist491 = []
dist492 = []
dist493 = []
dist494 = []
dist495 = []
dist496 = []
dist497 = []
dist498 = []
dist499 = []
dist500 = []
dist501 = []

#####substrate data#####
cnt=0
cnt1=0
cnt2=0
cnt3=0
for l in lines:
    if l[17:20] == 'UNL':
       l = l.strip('')
       #print(l)
       if l[21:22] == 'S':
           l = l.strip('')
           #print(l)
           if l[13:16] == 'N12':
                    l = l.strip('')
                    #print(l)
                    y15 = re.findall(r'[^\s]+', l)
                    Yvalue15.append(y15)

####protein data#####
for l in lines:
    if l[0:4] == 'ATOM':
        l = l.strip('')
        if l[21:22] == sys.argv[3]:
            l = l.strip('')
            #Arg
            if l[13:15] =='NE':
                l = l.strip('')
                x1 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue1.append(x1)
            if l[13:16] =='NH1' or l[13:16] == 'NH2':
                l = l.strip('')
                x2 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue2.append(x2)
            #Asn
            if l[13:16] =='ND2':
                l = l.strip('')
                x3 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue3.append(x3)
            #Asn
            if l[13:16] =='OD1':
                l = l.strip('')
                #x4 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                x4 = re.findall(r'[^\s]+', l)
                Xvalue4.append(x4)
            #Gln
            if l[13:16] =='NE2':
                l = l.strip('')
                x5 = re.findall(r'[^\s]+', l)
                #x5 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue5.append(x5)
            #Gln
            if  l[13:16] =='OE1':
                l = l.strip('')
                x6 = re.findall(r'[^\s]+', l)
                #x6 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue6.append(x6)
            #Glu
            if l[13:16] =='OE1' or l[13:16] =='OE2':
                l = l.strip('')
                #x1 = map(float,re.findall(r'[-+]?\d*\.\d+|\d+', l))
                x7 = re.findall(r'[^\s]+', l)
                Xvalue7.append(x7)
                cnt1 =cnt1 +1
            #Asp
            if l[13:16] =='OD1' or l[13:16] == 'OD2':
                l = l.strip('')
                x8 = re.findall(r'[^\s]+', l)
                #x2 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue8.append(x8)
            #HSD
            if l[13:16] =='ND1' or l[13:16] == 'NE2':
                l = l.strip('')
                x9 = re.findall(r'[^\s]+', l)
                #x7 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue9.append(x9)
            #Lys
            if l[13:15] =='NZ':
                l = l.strip('')
                x10 = re.findall(r'[^\s]+', l)
                #x8 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue10.append(x10)
            #Ser
            if l[13:15] =='OG':
                l = l.strip('')
                x11 = re.findall(r'[^\s]+', l)
                #x7 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue11.append(x11)
            #Thr
            if l[13:16] =='OG1':
                l = l.strip('')
                x12 = re.findall(r'[^\s]+', l)
                #x9 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue12.append(x12)
            #Trp
            if l[13:16] =='NE1':
                l = l.strip('')
                x13 = re.findall(r'[^\s]+', l)
                #x9 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue13.append(x13)
            #Tyr
            if l[13:15] =='OH':
                l = l.strip('')
                x14 = re.findall(r'[^\s]+', l)
                #x8 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Xvalue14.append(x14)
            #backbone O
            if l[13:14] == 'O':
                l = l.strip('')
                x15 = re.findall(r'[^\s]+', l)
                Xvalue15.append(x15)
    if l[0:4] == 'ATOM':
        l = l.strip('')
        if l[21:22] == sys.argv[4]:
            l = l.strip('')
            #Arg
            if l[13:15] =='NE':
                l = l.strip('')
                y1 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue1.append(y1)
            if l[13:16] =='NH1' or l[13:16] == 'NH2':
                l = l.strip('')
                y2 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue2.append(y2)
            #Asn
            if l[13:16] =='ND2':
                l = l.strip('')
                y3 = re.findall(r'[^\s]+', l)
                #x3 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue3.append(y3)
            #Asn
            if l[13:16] =='OD1':
                l = l.strip('')
                #x4 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                y4 = re.findall(r'[^\s]+', l)
                Yvalue4.append(y4)
            #Gln
            if l[13:16] =='NE2':
                l = l.strip('')
                y5 = re.findall(r'[^\s]+', l)
                #x5 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue5.append(y5)
            #Gln
            if  l[13:16] =='OE1':
                l = l.strip('')
                y6 = re.findall(r'[^\s]+', l)
                #x6 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue6.append(y6)
            #Glu
            if l[13:16] =='OE1' or l[13:16] =='OE2':
                l = l.strip('')
                #x1 = map(float,re.findall(r'[-+]?\d*\.\d+|\d+', l))
                y7 = re.findall(r'[^\s]+', l)
                Yvalue7.append(y7)
            #Asp
            if l[13:16] =='OD1' or l[13:16] == 'OD2':
                l = l.strip('')
                y8 = re.findall(r'[^\s]+', l)
                #x2 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue8.append(y8)
            #HSD
            if l[13:16] =='ND1' or l[13:16] == 'NE2':
                l = l.strip('')
                y9 = re.findall(r'[^\s]+', l)
                #x7 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue9.append(y9)
            #Lys
            if l[13:15] =='NZ':
                l = l.strip('')
                y10 = re.findall(r'[^\s]+', l)
                #x8 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue10.append(y10)
            #Ser
            if l[13:15] =='OG':
                l = l.strip('')
                y11 = re.findall(r'[^\s]+', l)
                #x7 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue11.append(y11)
            #Thr
            if l[13:16] =='OG1':
                l = l.strip('')
                y12 = re.findall(r'[^\s]+', l)
                #x9 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue12.append(y12)
            #Trp
            if l[13:16] =='NE1':
                l = l.strip('')
                y13 = re.findall(r'[^\s]+', l)
                #x9 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue13.append(y13)
            #Tyr
            if l[13:15] =='OH':
                l = l.strip('')
                y14 = re.findall(r'[^\s]+', l)
                #x8 = map(float, re.findall(r'[-+]?\d*\.\d+|\d+', l))
                Yvalue14.append(y14)
            #backbone N
            if l[13:14] == 'N':
                l = l.strip('')
                y16 = re.findall(r'[^\s]+', l)
                #y16 = re.findall(r'[A-Z]+', l)
                Yvalue16.append(y16)
            #backbone O
            if l[13:14] == 'O':
                l = l.strip('')
                y17 = re.findall(r'[^\s]+', l)
                #y17 = re.findall(r'[A-Z]+', l)
                #y17 = re.findall(r'[^\s]+', l)
                Yvalue17.append(y17)


####This starts the main loops that scans each array for  Xvalue or Yvalue 
###all scanned values are computed for sqrt(x^2 - y^2 -z^2)
###Then stored in new combined array

#ASNOD1-DONOR type
for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_1 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_1 = [float(Xvalue4[cnt4][5]), float(Yvalue1[cnt5][5]), dist_1]
                # print a
                dist1.append(a_1)

#dist_1 = np.linalg.norm(a-b)
for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_2 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_2 = [float(Xvalue4[cnt4][5]), float(Yvalue2[cnt5][5]), dist_2]
                # print a
                dist2.append(a_2)


for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_3 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_3 = [float(Xvalue4[cnt4][5]), float(Yvalue2[cnt5][5]), dist_3]
                # print a
                dist3.append(a_3)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_4 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_4 = [float(Xvalue4[cnt4][5]), float(Yvalue3[cnt5][5]), dist_4]
                # print a
                dist4.append(a_4)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_5 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_5 = [float(Xvalue4[cnt4][5]), float(Yvalue5[cnt5][5]), dist_5]
                # print a
                dist5.append(a_5)


for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_6 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_6 = [float(Xvalue4[cnt4][5]), float(Yvalue9[cnt5][5]), dist_6]
                # print a
                dist6.append(a_6)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_7 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_7 = [float(Xvalue4[cnt4][5]), float(Yvalue9[cnt5][5]), dist_7]
                # print a
                dist7.append(a_7)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_8 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_8 = [float(Xvalue4[cnt4][5]), float(Yvalue10[cnt5][5]), dist_8]
                # print a
                dist8.append(a_8)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_9 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_9 = [float(Xvalue4[cnt4][5]), float(Yvalue11[cnt5][5]), dist_9]
                # print a
                dist9.append(a_9)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_10 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_10 = [float(Xvalue4[cnt4][5]), float(Yvalue12[cnt5][5]), dist_10]
                # print a
                dist10.append(a_10)

for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_11 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_11 = [float(Xvalue4[cnt4][5]), float(Yvalue13[cnt5][5]), dist_11]
                # print a
                dist11.append(a_11)

#ASNOD1-DONOR type to chemical ligand (side-chain to side-chain) END of side chain to side chain
for cnt4 in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_12 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue4[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue4[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_12 = [float(Xvalue4[cnt4][5]), float(Yvalue14[cnt5][5]), dist_12]
                # print a
                dist12.append(a_12)



#############################################additional #######################################################
################################################################################################################

#ASNOD1-DONOR type to chemical ligand (side-chain to ligand)
for cnt4a in range(len(Xvalue4)):
    if Xvalue4[cnt4][3] == 'ASN' and Xvalue4[cnt4][4] == sys.argv[3] and Xvalue4[cnt4][2] == 'OD1':
        for cnt5a in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_13 = np.sqrt((float(Xvalue4[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue4[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue4[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_13 = [float(Xvalue4[cnt4][5]), float(Yvalue15[cnt5][5]), dist_13]
                # print a
                dist13.append(a_13)

#ASNOD1-backbone to chemical ligand (backbone to ligand) x
for cnt4b in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5b in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_14 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_14 = [float(Xvalue15[cnt4][5]), float(Yvalue15[cnt5][5]), dist_14]
                # print a
                dist14.append(a_14)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4c in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5c in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLN' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_15 = np.sqrt((float(Xvalue15[cnt4c][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_15 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_15]
                # print a
                dist15.append(a_15)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_16 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_16 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_16]
                # print a
                dist16.append(a_16)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_17 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_17 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_17]
                # print a
                dist17.append(a_17)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_18 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_18 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_18]
                # print a
                dist18.append(a_18)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_19 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_19 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_19]
                # print a
                dist19.append(a_19)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_20 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_20 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_20]
                # print a
                dist20.append(a_20)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_21 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_21 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_21]
                # print a
                dist21.append(a_21)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_22 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_22 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_22]
                # print a
                dist22.append(a_22)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_23 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_23 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_23]
                # print a
                dist23.append(a_23)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_24 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_24 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_24]
                # print a
                dist24.append(a_24)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_25 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_25 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_25]
                # print a
                dist25.append(a_25)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_26 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_26 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_26]
                # print a
                dist26.append(a_26)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_27 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_27 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_27]
                # print a
                dist27.append(a_27)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_28 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_28 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_28]
                # print a
                dist28.append(a_28)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_29 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_29 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_29]
                # print a
                dist29.append(a_29)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_30 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_30 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_30]
                # print a
                dist30.append(a_30)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_31 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_31 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_31]
                # print a
                dist31.append(a_31)

#ASNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_32 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_32 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_32]
                # print a
                dist32.append(a_32)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                dist_33 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_33 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_33]
                # print a
                dist33.append(a_33)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                dist_34 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_34 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_34]
                # print a
                dist34.append(a_34)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                dist_35 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_35 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_35]
                # print a
                dist35.append(a_35)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_36 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_36 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_36]
                # print a
                dist36.append(a_36)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_37 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_37 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_37]
                # print a
                dist37.append(a_37)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_38 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_38 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_38]
                # print a
                dist38.append(a_38)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_39 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_39 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_39]
                # print a
                dist39.append(a_39)

#ASNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == 'A' and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == 'A' and Yvalue10[cnt5][2] == 'NZ':
                dist_40 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_40 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_40]
                # print a
                dist40.append(a_40)

#ASNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_41 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_41 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_41]
                # print a
                dist41.append(a_41)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_42 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_42 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_42]
                # print a
                dist42.append(a_42)

#ASNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_43 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_43 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_43]
                # print a
                dist43.append(a_43)

#ASNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_44 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_44 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_44]
                # print a
                dist44.append(a_44)


###############################################################next##########################################
###############################################################################################################

#GLNOE1-DONOR type
for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_45 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_45 = [float(Xvalue6[cnt4][5]), float(Yvalue1[cnt5][5]), dist_45]
                # print a
                dist45.append(a_45)

#dist_1 = np.linalg.norm(a-b)
for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_46 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_46 = [float(Xvalue6[cnt4][5]), float(Yvalue2[cnt5][5]), dist_46]
                # print a
                dist46.append(a_46)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_47 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_47 = [float(Xvalue6[cnt4][5]), float(Yvalue2[cnt5][5]), dist_47]
                # print a
                dist47.append(a_47)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_48 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_48 = [float(Xvalue6[cnt4][5]), float(Yvalue3[cnt5][5]), dist_48]
                # print a
                dist48.append(a_48)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                dist_49 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_49 = [float(Xvalue6[cnt4][5]), float(Yvalue5[cnt5][5]), dist_49]
                # print a
                dist49.append(a_49)


for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                dist_50 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_50 = [float(Xvalue6[cnt4][5]), float(Yvalue9[cnt5][5]), dist_50]
                # print a
                dist50.append(a_50)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                dist_51 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_51 = [float(Xvalue6[cnt4][5]), float(Yvalue9[cnt5][5]), dist_51]
                # print a
                dist51.append(a_51)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                dist_52 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_52 = [float(Xvalue6[cnt4][5]), float(Yvalue10[cnt5][5]), dist_52]
                # print a
                dist52.append(a_52)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_53 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_53 = [float(Xvalue6[cnt4][5]), float(Yvalue11[cnt5][5]), dist_53]
                # print a
                dist53.append(a_53)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_54 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_54 = [float(Xvalue6[cnt4][5]), float(Yvalue12[cnt5][5]), dist_54]
                # print a
                dist54.append(a_54)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_55 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_55 = [float(Xvalue6[cnt4][5]), float(Yvalue13[cnt5][5]), dist_55]
                # print a
                dist55.append(a_55)

for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_56 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue6[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue6[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_56 = [float(Xvalue6[cnt4][5]), float(Yvalue14[cnt5][5]), dist_56]
                # print a
                dist56.append(a_56)



#################################### additional ##################################################################
##################################################################################################################


#GLNOD1-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue6)):
    if Xvalue6[cnt4][3] == 'GLN' and Xvalue6[cnt4][4] == sys.argv[3] and Xvalue6[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_57 = np.sqrt((float(Xvalue6[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue6[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue6[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_57 = [float(Xvalue6[cnt4][5]), float(Yvalue15[cnt5][5]), dist_57]
                # print a
                dist57.append(a_57)



#GLNOD1-backbone to chemical ligand (backbone to ligand) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_58 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_58 = [float(Xvalue15[cnt4][5]), float(Yvalue15[cnt5][5]), dist_58]
                # print a
                dist58.append(a_58)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLN' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_59 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_59 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_59]
                # print a
                dist59.append(a_59)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_60 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_60 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_60]
                # print a
                dist60.append(a_60)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_61 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_61 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_61]
                # print a
                dist61.append(a_61)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_62 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_62 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_62]
                # print a
                dist62.append(a_62)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_63 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_63 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_63]
                # print a
                dist63.append(a_63)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_64 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_64 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_64]
                # print a
                dist64.append(a_64)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_65 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_65 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_65]
                # print a
                dist65.append(a_65)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_66 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_66 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_66]
                # print a
                dist66.append(a_66)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_67 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_67 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_67]
                # print a
                dist67.append(a_67)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_68 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_68 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_68]
                # print a
                dist68.append(a_68)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_69 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_69 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_69]
                # print a
                dist69.append(a_69)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_70 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_70 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_70]
                # print a
                dist70.append(a_70)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_71 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_71 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_71]
                # print a
                dist71.append(a_71)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_72 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_72 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_72]
                # print a
                dist72.append(a_72)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_73 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_73 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_73]
                # print a
                dist73.append(a_73)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_74 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_74 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_74]
                # print a
                dist74.append(a_74)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_75 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_75 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_75]
                # print a
                dist75.append(a_75)

#GLNOD1-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_76 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_76 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_76]
                # print a
                dist76.append(a_76)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                dist_77 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_77 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_77]
                # print a
                dist77.append(a_77)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                dist_78 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_78 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_78]
                # print a
                dist78.append(a_78)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                dist_79 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_79 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_79]
                # print a
                dist79.append(a_79)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'GLN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_80 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_80 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_80]
                # print a
                dist80.append(a_80)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_81 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_81 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_81]
                # print a
                dist81.append(a_81)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_82 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_82 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_82]
                # print a
                dist82.append(a_82)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_83 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_83 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_83]
                # print a
                dist83.append(a_83)

#GLNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_84 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_84 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_84]
                # print a
                dist84.append(a_84)

#GLNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_85 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_85 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_85]
                # print a
                dist85.append(a_85)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_86 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_86 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_86]
                # print a
                dist86.append(a_86)

#GLNOD1-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_87 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_87 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_87]
                # print a
                dist87.append(a_87)

#GLNOD1-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLN' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_88 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_88 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_88]
                # print a
                dist88.append(a_88)




##############################################################################################################
##############################################################################################################

#GLUOE1-DONOR type
#GLUOE2-DONOR type
for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_89 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_89 = [float(Xvalue7[cnt4][5]), float(Yvalue1[cnt5][5]), dist_89]
                # print a
                dist89.append(a_89)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_90 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_90 = [float(Xvalue7[cnt4][5]), float(Yvalue1[cnt5][5]), dist_90]
                # print a
                dist90.append(a_90)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_91 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_91 = [float(Xvalue7[cnt4][5]), float(Yvalue2[cnt5][5]), dist_91]
                # print a
                dist91.append(a_91)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_92 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_92 = [float(Xvalue7[cnt4][5]), float(Yvalue2[cnt5][5]), dist_92]
                # print a
                dist92.append(a_92)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_93 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_93 = [float(Xvalue7[cnt4][5]), float(Yvalue2[cnt5][5]), dist_93]
                # print a
                dist93.append(a_93)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_94 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_94 = [float(Xvalue7[cnt4][5]), float(Yvalue2[cnt5][5]), dist_94]
                # print a
                dist94.append(a_94)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_95 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_95 = [float(Xvalue7[cnt4][5]), float(Yvalue3[cnt5][5]), dist_95]
                # print a
                dist95.append(a_95)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_96 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_96 = [float(Xvalue7[cnt4][5]), float(Yvalue3[cnt5][5]), dist_96]
                # print a
                dist96.append(a_96)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_97 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_97 = [float(Xvalue7[cnt4][5]), float(Yvalue5[cnt5][5]), dist_97]
                # print a
                dist97.append(a_97)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_98 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_98 = [float(Xvalue7[cnt4][5]), float(Yvalue5[cnt5][5]), dist_98]
                # print a
                dist98.append(a_98)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_99 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_99 = [float(Xvalue7[cnt4][5]), float(Yvalue9[cnt5][5]), dist_99]
                # print a
                dist99.append(a_99)

for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_100 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_100 = [float(Xvalue7[cnt4][5]), float(Yvalue9[cnt5][5]), dist_100]
                # print a
                dist100.append(a_100)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_101= np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_101 = [float(Xvalue7[cnt4][5]), float(Yvalue9[cnt5][5]), dist_101]
                # print a
                dist101.append(a_101)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_102 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_102 = [float(Xvalue7[cnt4][5]), float(Yvalue9[cnt5][5]), dist_102]
                # print a
                dist102.append(a_102)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_103 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_103 = [float(Xvalue7[cnt4][5]), float(Yvalue10[cnt5][5]), dist_103]
                # print a
                dist103.append(a_103)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_104 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_104 = [float(Xvalue7[cnt4][5]), float(Yvalue10[cnt5][5]), dist_104]
                # print a
                dist104.append(a_104)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_105 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_105 = [float(Xvalue7[cnt4][5]), float(Yvalue11[cnt5][5]), dist_105]
                # print a
                dist105.append(a_105)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_106 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_106 = [float(Xvalue7[cnt4][5]), float(Yvalue11[cnt5][5]), dist_106]
                # print a
                dist106.append(a_106)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_107 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_107 = [float(Xvalue7[cnt4][5]), float(Yvalue12[cnt5][5]), dist_107]
                # print a
                dist107.append(a_107)



for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_108 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_108 = [float(Xvalue7[cnt4][5]), float(Yvalue12[cnt5][5]), dist_108]
                # print a
                dist108.append(a_108)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_109 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_109 = [float(Xvalue7[cnt4][5]), float(Yvalue13[cnt5][5]), dist_109]
                # print a
                dist109.append(a_109)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_110 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_110 = [float(Xvalue7[cnt4][5]), float(Yvalue13[cnt5][5]), dist_110]
                # print a
                dist110.append(a_110)


for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_111 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_111 = [float(Xvalue7[cnt4][5]), float(Yvalue14[cnt5][5]), dist_111]
                # print a
                dist111.append(a_111)

for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_112 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue7[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue7[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_112 = [float(Xvalue7[cnt4][5]), float(Yvalue14[cnt5][5]), dist_112]
                # print a
                dist112.append(a_112)



################################################## additional ################################################
################################################################################################################


#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_113 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue7[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue7[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_113 = [float(Xvalue7[cnt4][5]), float(Yvalue15[cnt5][5]), dist_113]
                # print a
                dist113.append(a_113)


#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'GLU' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue16[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_114 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue7[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue7[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_114 = [float(Xvalue7[cnt4][5]), float(Yvalue15[cnt5][5]), dist_114]
                # print a
                dist114.append(a_114)


#-backbone to chemical ligand (backbone to ligand) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_115 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_115 = [float(Xvalue15[cnt4][5]), float(Yvalue15[cnt5][5]), dist_115]
                # print a
                dist115.append(a_115)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_116 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_116 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_116]
                # print a
                dist116.append(a_116)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'GLU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_117 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_117 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_117]
                # print a
                dist117.append(a_117)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_118 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_118 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_118]
                # print a
                dist118.append(a_118)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_119 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_119 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_119]
                # print a
                dist119.append(a_119)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_120 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_120 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_120]
                # print a
                dist120.append(a_120)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_121 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_121 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_121]
                # print a
                dist121.append(a_121)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_122 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_122 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_122]
                # print a
                dist122.append(a_122)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_123 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_123 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_123]
                # print a
                dist23.append(a_123)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_124 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_124 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_124]
                # print a
                dist124.append(a_124)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_125 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_125 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_125]
                # print a
                dist125.append(a_125)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_126 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_126 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_126]
                # print a
                dist126.append(a_126)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_127 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_127 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_127]
                # print a
                dist127.append(a_127)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_128 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_128 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_128]
                # print a
                dist128.append(a_128)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_129 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_129 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_129]
                # print a
                dist129.append(a_129)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_130 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_130 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_130]
                # print a
                dist130.append(a_130)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_131 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_131 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_131]
                # print a
                dist131.append(a_131)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_132 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_132 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_132]
                # print a
                dist132.append(a_132)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_133 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_133 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_133]
                # print a
                dist133.append(a_133)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                dist_134 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_134 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_134]
                # print a
                dist134.append(a_134)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                dist_135 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_135 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_135]
                # print a
                dist135.append(a_135)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                dist_136 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_136 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_136]
                # print a
                dist136.append(a_136)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'GLU' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_137 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_137 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_137]
                # print a
                dist137.append(a_137)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLU' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_138 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_138 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_138]
                # print a
                dist138.append(a_138)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_139 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_139 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_139]
                # print a
                dist139.append(a_139)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_140 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_140 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_140]
                # print a
                dist140.append(a_140)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_141 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_141 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_141]
                # print a
                dist141.append(a_141)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_142 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_142 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_142]
                # print a
                dist142.append(a_142)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_142a = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_142a = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_142a]
                # print a
                dist142a.append(a_142a)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_143 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_143 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_143]
                # print a
                dist143.append(a_143)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_144 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_144 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_144]
                # print a
                dist144.append(a_144)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'GLU' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_145 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_145 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_145]
                # print a
                dist145.append(a_145)

##############################################################################################################

#ASPOD1-DONOR type
#ASPOD2-DONOR type
for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_146 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_146 = [float(Xvalue8[cnt4][5]), float(Yvalue1[cnt5][5]), dist_146]
                # print a
                dist146.append(a_146)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_147 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_147 = [float(Xvalue8[cnt4][5]), float(Yvalue1[cnt5][5]), dist_147]
                # print a
                dist147.append(a_147)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_148 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_148 = [float(Xvalue8[cnt4][5]), float(Yvalue2[cnt5][5]), dist_148]
                # print a
                dist148.append(a_148)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_149 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_149 = [float(Xvalue8[cnt4][5]), float(Yvalue2[cnt5][5]), dist_149]
                # print a
                dist149.append(a_149)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_150 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_150 = [float(Xvalue8[cnt4][5]), float(Yvalue2[cnt5][5]), dist_150]
                # print a
                dist150.append(a_150)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_151 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_151 = [float(Xvalue8[cnt4][5]), float(Yvalue2[cnt5][5]), dist_151]
                # print a
                dist151.append(a_151)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_152 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_152 = [float(Xvalue8[cnt4][5]), float(Yvalue3[cnt5][5]), dist_152]
                # print a
                dist152.append(a_152)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_153 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_153 = [float(Xvalue8[cnt4][5]), float(Yvalue3[cnt5][5]), dist_153]
                # print a
                dist153.append(a_153)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_154 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_154 = [float(Xvalue8[cnt4][5]), float(Yvalue5[cnt5][5]), dist_154]
                # print a
                dist154.append(a_154)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_155 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_155 = [float(Xvalue8[cnt4][5]), float(Yvalue5[cnt5][5]), dist_155]
                # print a
                dist155.append(a_155)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_156 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_156 = [float(Xvalue8[cnt4][5]), float(Yvalue9[cnt5][5]), dist_156]
                # print a
                dist156.append(a_156)

for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_157 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_157 = [float(Xvalue8[cnt4][5]), float(Yvalue9[cnt5][5]), dist_157]
                # print a
                dist157.append(a_157)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_158= np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_158 = [float(Xvalue8[cnt4][5]), float(Yvalue9[cnt5][5]), dist_158]
                # print a
                dist158.append(a_158)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_159 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_159 = [float(Xvalue8[cnt4][5]), float(Yvalue9[cnt5][5]), dist_159]
                # print a
                dist159.append(a_159)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_160 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_160 = [float(Xvalue8[cnt4][5]), float(Yvalue10[cnt5][5]), dist_160]
                # print a
                dist160.append(a_160)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_161 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_161 = [float(Xvalue8[cnt4][5]), float(Yvalue10[cnt5][5]), dist_161]
                # print a
                dist161.append(a_161)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_162 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_162 = [float(Xvalue8[cnt4][5]), float(Yvalue11[cnt5][5]), dist_162]
                # print a
                dist162.append(a_162)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_163 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_163 = [float(Xvalue8[cnt4][5]), float(Yvalue11[cnt5][5]), dist_163]
                # print a
                dist163.append(a_163)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_164 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_164= [float(Xvalue8[cnt4][5]), float(Yvalue12[cnt5][5]), dist_164]
                # print a
                dist164.append(a_164)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_165 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_165 = [float(Xvalue8[cnt4][5]), float(Yvalue12[cnt5][5]), dist_165]
                # print a
                dist165.append(a_165)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_166 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_166 = [float(Xvalue8[cnt4][5]), float(Yvalue13[cnt5][5]), dist_166]
                # print a
                dist166.append(a_166)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_167 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_167 = [float(Xvalue8[cnt4][5]), float(Yvalue13[cnt5][5]), dist_167]
                # print a
                dist167.append(a_167)



for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_168 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_168 = [float(Xvalue8[cnt4][5]), float(Yvalue14[cnt5][5]), dist_168]
                # print a
                dist168.append(a_168)


for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_169 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue8[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue8[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_169 = [float(Xvalue8[cnt4][5]), float(Yvalue14[cnt5][5]), dist_169]
                # print a
                dist169.append(a_169)


################################ADDITIONAL#############################################################
#########################################################################################################

#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue8)):
    if Xvalue8[cnt4][3] == 'ASP' and Xvalue8[cnt4][4] == sys.argv[3] and Xvalue8[cnt4][2] == 'OD2':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_170 = np.sqrt((float(Xvalue8[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue8[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue8[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_170 = [float(Xvalue8[cnt4][5]), float(Yvalue15[cnt5][5]), dist_170]
                # print a
                dist170.append(a_170)


#-backbone to chemical ligand (backbone to ligand) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_171 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_171 = [float(Xvalue15[cnt4][5]), float(Yvalue15[cnt5][5]), dist_171]
                # print a
                dist171.append(a_171)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_172 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_172 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_172]
                # print a
                dist172.append(a_172)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_173 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_173 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_173]
                # print a
                dist173.append(a_173)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ASP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_174 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_174 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_174]
                # print a
                dist174.append(a_174)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_175 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_175 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_175]
                # print a
                dist175.append(a_175)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_176 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_176 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_176]
                # print a
                dist176.append(a_176)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_177 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_177 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_177]
                # print a
                dist177.append(a_177)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_178 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_178 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_178]
                # print a
                dist178.append(a_178)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_179 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_179 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_179]
                # print a
                dist179.append(a_179)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_180 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_180 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_180]
                # print a
                dist180.append(a_180)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_181 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_181 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_181]
                # print a
                dist181.append(a_181)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_182 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_182 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_182]
                # print a
                dist182.append(a_182)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_183 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_183 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_183]
                # print a
                dist183.append(a_183)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_184 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_184 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_184]
                # print a
                dist184.append(a_184)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_185 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_185 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_185]
                # print a
                dist185.append(a_185)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_186 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_186 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_186]
                # print a
                dist186.append(a_186)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_187 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_187 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_187]
                # print a
                dist187.append(a_187)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_188 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_188 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_188]
                # print a
                dist188.append(a_188)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_189 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_189 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_189]
                # print a
                dist189.append(a_189)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[4] and Yvalue1[cnt5][2] == 'NE':
                dist_190 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_190 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_190]
                # print a
                dist190.append(a_190)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH1':
                dist_191 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_191 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_191]
                # print a
                dist191.append(a_191)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[4] and Yvalue2[cnt5][2] == 'NH2':
                dist_192 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_192 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_192]
                # print a
                dist192.append(a_192)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASP' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_193 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_193 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_193]
                # print a
                dist193.append(a_193)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'ASP' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_194 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_194 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_194]
                # print a
                dist194.append(a_194)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_195 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_195 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_195]
                # print a
                dist195.append(a_195)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_196 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_196 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_196]
                # print a
                dist196.append(a_196)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_197 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_197 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_197]
                # print a
                dist197.append(a_197)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_198 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_198 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_198]
                # print a
                dist198.append(a_198)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_198a = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_198a = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_198a]
                # print a
                dist198a.append(a_198a)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_199 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_199 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_199]
                # print a
                dist199.append(a_199)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_200 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_200 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_200]
                # print a
                dist200.append(a_200)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'ASP' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_201 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_201 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_201]
                # print a
                dist201.append(a_201)






############################################################################################################


#HISND1-DONOR type
#HISNE2-DONOR type
for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_202 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_202 = [float(Xvalue9[cnt4][5]), float(Yvalue1[cnt5][5]), dist_202]
                # print a
                dist202.append(a_202)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_203 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_203 = [float(Xvalue9[cnt4][5]), float(Yvalue1[cnt5][5]), dist_203]
                # print a
                dist203.append(a_203)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_204 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_204 = [float(Xvalue9[cnt4][5]), float(Yvalue2[cnt5][5]), dist_204]
                # print a
                dist204.append(a_204)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_205 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_205 = [float(Xvalue9[cnt4][5]), float(Yvalue2[cnt5][5]), dist_205]
                # print a
                dist205.append(a_205)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_206 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_206 = [float(Xvalue9[cnt4][5]), float(Yvalue2[cnt5][5]), dist_206]
                # print a
                dist206.append(a_206)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_207 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_207 = [float(Xvalue9[cnt4][5]), float(Yvalue2[cnt5][5]), dist_207]
                # print a
                dist207.append(a_207)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_208 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_208 = [float(Xvalue9[cnt4][5]), float(Yvalue3[cnt5][5]), dist_208]
                # print a
                dist208.append(a_208)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_209 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_209 = [float(Xvalue9[cnt4][5]), float(Yvalue3[cnt5][5]), dist_209]
                # print a
                dist209.append(a_209)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_210 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_210 = [float(Xvalue9[cnt4][5]), float(Yvalue5[cnt5][5]), dist_210]
                # print a
                dist210.append(a_210)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_211 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_211 = [float(Xvalue9[cnt4][5]), float(Yvalue5[cnt5][5]), dist_211]
                # print a
                dist211.append(a_211)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_212 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_212 = [float(Xvalue9[cnt4][5]), float(Yvalue9[cnt5][5]), dist_212]
                # print a
                dist212.append(a_212)

for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_213 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_213 = [float(Xvalue9[cnt4][5]), float(Yvalue9[cnt5][5]), dist_213]
                # print a
                dist213.append(a_213)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_214= np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_214 = [float(Xvalue9[cnt4][5]), float(Yvalue9[cnt5][5]), dist_214]
                # print a
                dist214.append(a_214)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_215 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_215 = [float(Xvalue9[cnt4][5]), float(Yvalue9[cnt5][5]), dist_215]
                # print a
                dist215.append(a_215)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_216 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_216 = [float(Xvalue9[cnt4][5]), float(Yvalue10[cnt5][5]), dist_216]
                # print a
                dist216.append(a_216)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_217 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_217 = [float(Xvalue9[cnt4][5]), float(Yvalue10[cnt5][5]), dist_217]
                # print a
                dist217.append(a_217)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_218 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_218 = [float(Xvalue9[cnt4][5]), float(Yvalue11[cnt5][5]), dist_18]
                # print a
                dist218.append(a_218)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_219 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_219 = [float(Xvalue9[cnt4][5]), float(Yvalue11[cnt5][5]), dist_219]
                # print a
                dist219.append(a_219)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_220 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_220 = [float(Xvalue9[cnt4][5]), float(Yvalue12[cnt5][5]), dist_220]
                # print a
                dist220.append(a_220)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_221 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_221 = [float(Xvalue9[cnt4][5]), float(Yvalue12[cnt5][5]), dist_221]
                # print a
                dist221.append(a_221)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_222 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_222 = [float(Xvalue9[cnt4][5]), float(Yvalue13[cnt5][5]), dist_222]
                # print a
                dist222.append(a_222)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_223 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_223 = [float(Xvalue9[cnt4][5]), float(Yvalue13[cnt5][5]), dist_223]
                # print a
                dist223.append(a_223)



for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_224 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_224 = [float(Xvalue9[cnt4][5]), float(Yvalue14[cnt5][5]), dist_224]
                # print a
                dist224.append(a_224)


for cnt4 in range(len(Xvalue9)):
    if Xvalue9[cnt4][3] == 'HSD' and Xvalue9[cnt4][4] == sys.argv[3] and Xvalue9[cnt4][2] == 'NE2':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_225 = np.sqrt((float(Xvalue9[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue9[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue9[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_225 = [float(Xvalue9[cnt4][5]), float(Yvalue14[cnt5][5]), dist_225]
                # print a
                dist225.append(a_225)

######################################################################################################################
##################################################ADDITIONAL##########################################################

#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue7)):
    if Xvalue7[cnt4][3] == 'HSD' and Xvalue7[cnt4][4] == sys.argv[3] and Xvalue7[cnt4][2] == 'ND1':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_226 = np.sqrt((float(Xvalue7[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue7[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue7[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_226 = [float(Xvalue8[cnt4][5]), float(Yvalue15[cnt5][5]), dist_226]
                # print a
                dist226.append(a_226)


#-backbone to chemical ligand (backbone to ligand) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_227 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_227 = [float(Xvalue15[cnt4][5]), float(Yvalue15[cnt5][5]), dist_227]
                # print a
                dist227.append(a_227)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_228 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_228 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_228]
                # print a
                dist228.append(a_228)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_229 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_229 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_229]
                # print a
                dist229.append(a_229)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_230 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_230 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_230]
                # print a
                dist230.append(a_230)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_231 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_231 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_231]
                # print a
                dist231.append(a_231)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_232 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_232 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_232]
                # print a
                dist232.append(a_232)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_233 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_233 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_233]
                # print a
                dist233.append(a_233)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_234 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_234 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_234]
                # print a
                dist234.append(a_234)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_235 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_235 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_235]
                # print a
                dist235.append(a_235)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_236 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_236 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_236]
                # print a
                dist236.append(a_236)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_237 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_237 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_237]
                # print a
                dist237.append(a_237)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_238 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_238 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_238]
                # print a
                dist238.append(a_238)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'HSD' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_239 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_239 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_239]
                # print a
                dist239.append(a_239)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_240 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_240 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_240]
                # print a
                dist240.append(a_240)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_241 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_241 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_241]
                # print a
                dist241.append(a_241)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_242 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_242 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_242]
                # print a
                dist242.append(a_242)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_243 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_243 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_243]
                # print a
                dist243.append(a_243)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == 'B' and Yvalue1[cnt5][2] == 'NE':
                dist_244 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_244 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_244]
                # print a
                dist244.append(a_244)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH1':
                dist_245 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_245 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_245]
                # print a
                dist245.append(a_245)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH2':
                dist_246 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_246 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_246]
                # print a
                dist246.append(a_246)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'HSD' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_247 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_247 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_247]
                # print a
                dist247.append(a_247)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'HSD' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_248 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_248 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_248]
                # print a
                dist248.append(a_248)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_249 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_249 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_249]
                # print a
                dist249.append(a_249)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_250 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_250 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_250]
                # print a
                dist250.append(a_250)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_251 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_251 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_251]
                # print a
                dist251.append(a_251)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_252 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_252 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_252]
                # print a
                dist252.append(a_252)



#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_253 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_253 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_253]
                # print a
                dist253.append(a_253)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_254 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_254 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_254]
                # print a
                dist254.append(a_254)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'HSD' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_255 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_255 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_255]
                # print a
                dist255.append(a_255)




##############################################################################################################
#SEROG-DONOR type
for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_256 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_256 = [float(Xvalue11[cnt4][5]), float(Yvalue1[cnt5][5]), dist_256]
                # print a
                dist256.append(a_256)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_257 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_257 = [float(Xvalue11[cnt4][5]), float(Yvalue2[cnt5][5]), dist_257]
                # print a
                dist257.append(a_257)

for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_258 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_258 = [float(Xvalue11[cnt4][5]), float(Yvalue2[cnt5][5]), dist_258]
                # print a
                dist258.append(a_258)

for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_259 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_259 = [float(Xvalue11[cnt4][5]), float(Yvalue3[cnt5][5]), dist_259]
                # print a
                dist259.append(a_259)

for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_260 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_260 = [float(Xvalue11[cnt4][5]), float(Yvalue5[cnt5][5]), dist_260]
                # print a
                dist260.append(a_260)

for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_261 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_261 = [float(Xvalue11[cnt4][5]), float(Yvalue9[cnt5][5]), dist_261]
                # print a
                dist261.append(a_261)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_262= np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_262 = [float(Xvalue11[cnt4][5]), float(Yvalue9[cnt5][5]), dist_262]
                # print a
                dist262.append(a_262)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_263 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_263 = [float(Xvalue11[cnt4][5]), float(Yvalue10[cnt5][5]), dist_263]
                # print a
                dist263.append(a_263)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_264 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_264 = [float(Xvalue11[cnt4][5]), float(Yvalue11[cnt5][5]), dist_264]
                # print a
                dist264.append(a_264)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_265 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_265 = [float(Xvalue11[cnt4][5]), float(Yvalue12[cnt5][5]), dist_265]
                # print a
                dist265.append(a_265)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_266 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_266 = [float(Xvalue11[cnt4][5]), float(Yvalue13[cnt5][5]), dist_266]
                # print a
                dist266.append(a_266)


for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_267 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue11[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue11[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_267 = [float(Xvalue11[cnt4][5]), float(Yvalue14[cnt5][5]), dist_267]
                # print a
                dist267.append(a_267)


######################ADDITIONAL###############################################################################
################################################################################################################
#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue11)):
    if Xvalue11[cnt4][3] == 'SER' and Xvalue11[cnt4][4] == sys.argv[3] and Xvalue11[cnt4][2] == 'OG':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_268 = np.sqrt((float(Xvalue11[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue11[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue11[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_268 = [float(Xvalue11[cnt4][5]), float(Yvalue15[cnt5][5]), dist_268]
                # print a
                dist268.append(a_268)



#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_269 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_269 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_269]
                # print a
                dist269.append(a_269)



#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_270 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_270 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_270]
                # print a
                dist270.append(a_270)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_271 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_271 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_271]
                # print a
                dist271.append(a_271)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_272 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_272 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_272]
                # print a
                dist272.append(a_272)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_273 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_273 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_273]
                # print a
                dist273.append(a_273)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_274 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_274 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_274]
                # print a
                dist274.append(a_274)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_275 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_275 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_275]
                # print a
                dist275.append(a_275)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_276 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_276 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_276]
                # print a
                dist276.append(a_276)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_277 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_277 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_277]
                # print a
                dist277.append(a_277)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_278 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_278 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_278]
                # print a
                dist278.append(a_278)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_279 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_279 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_279]
                # print a
                dist279.append(a_279)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'SER' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_280 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_280 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_280]
                # print a
                dist280.append(a_280)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_281 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_281 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_281]
                # print a
                dist281.append(a_281)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_282 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_282 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_282]
                # print a
                dist282.append(a_282)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_283 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_283 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_283]
                # print a
                dist283.append(a_283)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_284 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_284 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_284]
                # print a
                dist284.append(a_284)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == 'B' and Yvalue1[cnt5][2] == 'NE':
                dist_285 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_285 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_285]
                # print a
                dist285.append(a_285)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH1':
                dist_286 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_286 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_286]
                # print a
                dist286.append(a_286)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH2':
                dist_287 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_287 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_287]
                # print a
                dist287.append(a_287)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_288 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_288 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_288]
                # print a
                dist288.append(a_288)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_289 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_289 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_289]
                # print a
                dist289.append(a_289)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HIS' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_290 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_290 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_290]
                # print a
                dist290.append(a_290)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HIS' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_291 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_291 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_291]
                # print a
                dist291.append(a_291)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_292 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_292 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_292]
                # print a
                dist292.append(a_292)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_293 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_293 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_293]
                # print a
                dist293.append(a_293)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_293a = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_293a = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_293a]
                # print a
                dist293a.append(a_293a)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_294 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_294 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_294]
                # print a
                dist294.append(a_294)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_295 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_295 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_295]
                # print a
                dist295.append(a_295)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'SER' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_296 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_296 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_296]
                # print a
                dist296.append(a_296)

###############################################################################################################


#THROG11-DONOR type
for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_297 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_297 = [float(Xvalue12[cnt4][5]), float(Yvalue1[cnt5][5]), dist_297]
                # print a
                dist297.append(a_297)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_298 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_298 = [float(Xvalue12[cnt4][5]), float(Yvalue2[cnt5][5]), dist_298]
                # print a
                dist298.append(a_298)

for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_299 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_299 = [float(Xvalue12[cnt4][5]), float(Yvalue2[cnt5][5]), dist_299]
                # print a
                dist299.append(a_299)

for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_300 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_300 = [float(Xvalue12[cnt4][5]), float(Yvalue3[cnt5][5]), dist_300]
                # print a
                dist300.append(a_300)

for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_301 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_301 = [float(Xvalue12[cnt4][5]), float(Yvalue5[cnt5][5]), dist_301]
                # print a
                dist301.append(a_301)

for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_302 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_302 = [float(Xvalue12[cnt4][5]), float(Yvalue9[cnt5][5]), dist_302]
                # print a
                dist302.append(a_302)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_303= np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_303 = [float(Xvalue12[cnt4][5]), float(Yvalue9[cnt5][5]), dist_303]
                # print a
                dist303.append(a_303)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_304 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_304 = [float(Xvalue12[cnt4][5]), float(Yvalue10[cnt5][5]), dist_304]
                # print a
                dist304.append(a_304)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_305 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_305 = [float(Xvalue12[cnt4][5]), float(Yvalue11[cnt5][5]), dist_305]
                # print a
                dist305.append(a_305)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_306 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_306 = [float(Xvalue12[cnt4][5]), float(Yvalue12[cnt5][5]), dist_306]
                # print a
                dist306.append(a_306)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_307 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_307 = [float(Xvalue12[cnt4][5]), float(Yvalue13[cnt5][5]), dist_307]
                # print a
                dist307.append(a_307)


for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_308 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue12[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue12[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_308 = [float(Xvalue12[cnt4][5]), float(Yvalue14[cnt5][5]), dist_308]
                # print a
                dist308.append(a_308)
##########################################################################################

#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue12)):
    if Xvalue12[cnt4][3] == 'THR' and Xvalue12[cnt4][4] == sys.argv[3] and Xvalue12[cnt4][2] == 'OG1':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_309 = np.sqrt((float(Xvalue12[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue12[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue12[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_309 = [float(Xvalue12[cnt4][5]), float(Yvalue15[cnt5][5]), dist_309]
                # print a
                dist309.append(a_309)



#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_310 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_310 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_310]
                # print a
                dist310.append(a_310)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_311 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_311 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_311]
                # print a
                dist311.append(a_311)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_312 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_312 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_312]
                # print a
                dist312.append(a_312)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_313 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_313 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_313]
                # print a
                dist313.append(a_313)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_314 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_314 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_314]
                # print a
                dist314.append(a_314)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_315 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_315 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_315]
                # print a
                dist315.append(a_315)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_316 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_316 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_316]
                # print a
                dist316.append(a_316)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_317 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_317 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_317]
                # print a
                dist317.append(a_317)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_318 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_318 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_318]
                # print a
                dist318.append(a_318)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_319 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_319 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_319]
                # print a
                dist319.append(a_319)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_320 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_320 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_320]
                # print a
                dist320.append(a_320)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_321 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_321 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_321]
                # print a
                dist321.append(a_321)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_322 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_322 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_322]
                # print a
                dist322.append(a_322)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_323 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_323 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_323]
                # print a
                dist323.append(a_323)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_324 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_324 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_324]
                # print a
                dist324.append(a_324)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_325 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_325 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_325]
                # print a
                dist325.append(a_325)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'THR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_326 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_326 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_326]
                # print a
                dist326.append(a_326)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == 'B' and Yvalue1[cnt5][2] == 'NE':
                dist_327 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_327 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_327]
                # print a
                dist327.append(a_327)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH1':
                dist_328 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_328 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_328]
                # print a
                dist328.append(a_328)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH2':
                dist_329 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_329 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_329]
                # print a
                dist329.append(a_329)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'THR' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_330 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_330 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_330]
                # print a
                dist330.append(a_330)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'THR' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_331 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_331 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_331]
                # print a
                dist331.append(a_331)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'THR' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_332 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_332 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_332]
                # print a
                dist332.append(a_332)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'THR' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_333 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_333 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_333]
                # print a
                dist333.append(a_333)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_334 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_334 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_334]
                # print a
                dist334.append(a_334)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE1':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'THR' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_335 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_335 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_335]
                # print a
                dist335.append(a_335)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'OE2':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'THR' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_336 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_336 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_336]
                # print a
                dist336.append(a_336)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_337 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_337 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_337]
                # print a
                dist337.append(a_337)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_338 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_338 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_338]
                # print a
                dist338.append(a_338)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'THR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_339 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_339 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_339]
                # print a
                dist339.append(a_339)

                
############################################################################################
#TYROH-DONOR type
for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == sys.argv[3] and Yvalue1[cnt5][2] == 'NE':
                #print cnt4,cnt5
                dist_340 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_340 = [float(Xvalue14[cnt4][5]), float(Yvalue1[cnt5][5]), dist_340]
                # print a
                dist340.append(a_340)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH1':
                #print cnt4,cnt5
                dist_341 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_341 = [float(Xvalue14[cnt4][5]), float(Yvalue2[cnt5][5]), dist_341]
                # print a
                dist341.append(a_341)

for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == sys.argv[3] and Yvalue2[cnt5][2] == 'NH2':
                #print cnt4,cnt5
                dist_342 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_342 = [float(Xvalue14[cnt4][5]), float(Yvalue2[cnt5][5]), dist_342]
                # print a
                dist342.append(a_342)

for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'ASN' and Yvalue3[cnt5][4] == sys.argv[3] and Yvalue3[cnt5][2] == 'ND2':
                #print cnt4,cnt5
                dist_343 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_343 = [float(Xvalue14[cnt4][5]), float(Yvalue3[cnt5][5]), dist_343]
                # print a
                dist343.append(a_343)

for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'GLN' and Yvalue5[cnt5][4] == sys.argv[3] and Yvalue5[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_344 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_344 = [float(Xvalue14[cnt4][5]), float(Yvalue5[cnt5][5]), dist_344]
                # print a
                dist344.append(a_344)

for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'ND1':
                #print cnt4,cnt5
                dist_345 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_345 = [float(Xvalue14[cnt4][5]), float(Yvalue9[cnt5][5]), dist_345]
                # print a
                dist345.append(a_345)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'HSD' and Yvalue9[cnt5][4] == sys.argv[3] and Yvalue9[cnt5][2] == 'NE2':
                #print cnt4,cnt5
                dist_346= np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_346 = [float(Xvalue14[cnt4][5]), float(Yvalue9[cnt5][5]), dist_346]
                # print a
                dist346.append(a_346)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[3] and Yvalue10[cnt5][2] == 'NZ':
                #print cnt4,cnt5
                dist_347 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_347 = [float(Xvalue14[cnt4][5]), float(Yvalue10[cnt5][5]), dist_347]
                # print a
                dist347.append(a_347)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[3] and Yvalue11[cnt5][2] == 'OG':
                #print cnt4,cnt5
                dist_348 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_348 = [float(Xvalue14[cnt4][5]), float(Yvalue11[cnt5][5]), dist_348]
                # print a
                dist348.append(a_348)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[3] and Yvalue12[cnt5][2] == 'OG1':
                #print cnt4,cnt5
                dist_349 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_349 = [float(Xvalue14[cnt4][5]), float(Yvalue12[cnt5][5]), dist_349]
                # print a
                dist349.append(a_349)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[3] and Yvalue13[cnt5][2] == 'NE1':
                #print cnt4,cnt5
                dist_350 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_350 = [float(Xvalue14[cnt4][5]), float(Yvalue13[cnt5][5]), dist_350]
                # print a
                dist350.append(a_350)


for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[3] and Yvalue14[cnt5][2] == 'OH':
                #print cnt4,cnt5
                dist_351 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                                float(Xvalue14[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                               float(Xvalue14[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_351 = [float(Xvalue14[cnt4][5]), float(Yvalue14[cnt5][5]), dist_351]
                # print a
                dist351.append(a_351)

####################################
#-DONOR type to chemical ligand (side-chain to ligand)
for cnt4 in range(len(Xvalue14)):
    if Xvalue14[cnt4][3] == 'TYR' and Xvalue14[cnt4][4] == sys.argv[3] and Xvalue14[cnt4][2] == 'OH':
        for cnt5 in range(len(Yvalue15)):
            if Yvalue15[cnt5][3] == 'UNL' and Yvalue15[cnt5][4] == 'S' and Yvalue15[cnt5][2] == 'N12':
                dist_352 = np.sqrt((float(Xvalue14[cnt4][6]) - float(Yvalue15[cnt5][6])) ** 2 + (
                        float(Xvalue14[cnt4][7]) - float(Yvalue15[cnt5][7])) ** 2 + (
                                         float(Xvalue14[cnt4][8]) - float(Yvalue15[cnt5][8])) ** 2)
                a_352 = [float(Xvalue14[cnt4][5]), float(Yvalue15[cnt5][5]), dist_352]
                # print a
                dist352.append(a_352)



#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_353= np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_353 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_353]
                # print a
                dist353.append(a_353)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_354 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_354 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_354]
                # print a
                dist354.append(a_354)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_355 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_355 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_355]
                # print a
                dist355.append(a_355)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'VAL' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_356 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_356 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_356]
                # print a
                dist356.append(a_356)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ILE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_357 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_357 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_357]
                # print a
                dist357.append(a_357)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LEU' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_358 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_358 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_358]
                # print a
                dist358.append(a_358)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'MET' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_359 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_359 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_359]
                # print a
                dist359.append(a_359)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_360 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_360 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_360]
                # print a
                dist360.append(a_360)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'CYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_361 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_361 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_361]
                # print a
                dist361.append(a_361)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TRP' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_362 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_362 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_362]
                # print a
                dist362.append(a_362)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_363 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_363 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_363]
                # print a
                dist363.append(a_363)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PRO' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_364 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_364 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_364]
                # print a
                dist364.append(a_364)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ALA' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_365 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_365 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_365]
                # print a
                dist365.append(a_365)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_366 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_366 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_366]
                # print a
                dist366.append(a_366)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'ARG' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_367 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_367 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_367]
                # print a
                dist367.append(a_367)

#backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'LYS' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_368 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_368 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_368]
                # print a
                dist368.append(a_368)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'PHE' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_369 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_369 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_369]
                # print a
                dist369.append(a_369)

#-backbone to chemical ligand (backbone to backbone) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue16)):
            if Yvalue16[cnt5][3] == 'TYR' and Yvalue16[cnt5][4] == sys.argv[4] and Yvalue16[cnt5][2] == 'N':
                dist_370 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue16[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue16[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue16[cnt5][8])) ** 2)
                a_370 = [float(Xvalue15[cnt4][5]), float(Yvalue16[cnt5][5]), dist_370]
                # print a
                dist370.append(a_370)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue1)):
            if Yvalue1[cnt5][3] == 'ARG' and Yvalue1[cnt5][4] == 'B' and Yvalue1[cnt5][2] == 'NE':
                dist_371 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue1[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue1[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue1[cnt5][8])) ** 2)
                a_371 = [float(Xvalue15[cnt4][5]), float(Yvalue1[cnt5][5]), dist_371]
                # print a
                dist371.append(a_371)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH1':
                dist_372 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_372 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_372]
                # print a
                dist372.append(a_372)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue2)):
            if Yvalue2[cnt5][3] == 'ARG' and Yvalue2[cnt5][4] == 'B' and Yvalue2[cnt5][2] == 'NH2':
                dist_373 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue2[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue2[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue2[cnt5][8])) ** 2)
                a_373 = [float(Xvalue15[cnt4][5]), float(Yvalue2[cnt5][5]), dist_373]
                # print a
                dist373.append(a_373)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue3)):
            if Yvalue3[cnt5][3] == 'TYR' and Yvalue3[cnt5][4] == sys.argv[4] and Yvalue3[cnt5][2] == 'ND2':
                dist_374 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue3[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue3[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue3[cnt5][8])) ** 2)
                a_374 = [float(Xvalue15[cnt4][5]), float(Yvalue3[cnt5][5]), dist_374]
                # print a
                dist374.append(a_374)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue5)):
            if Yvalue5[cnt5][3] == 'TYR' and Yvalue5[cnt5][4] == sys.argv[4] and Yvalue5[cnt5][2] == 'NE2':
                dist_375 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue5[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue5[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue5[cnt5][8])) ** 2)
                a_375 = [float(Xvalue15[cnt4][5]), float(Yvalue5[cnt5][5]), dist_375]
                # print a
                dist375.append(a_375)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'TYR' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'ND1':
                dist_376 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_376 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_376]
                # print a
                dist376.append(a_376)

#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue9)):
            if Yvalue9[cnt5][3] == 'TYR' and Yvalue9[cnt5][4] == sys.argv[4] and Yvalue9[cnt5][2] == 'NE2':
                dist_377 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue9[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue9[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue9[cnt5][8])) ** 2)
                a_377 = [float(Xvalue15[cnt4][5]), float(Yvalue9[cnt5][5]), dist_377]
                # print a
                dist377.append(a_377)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue10)):
            if Yvalue10[cnt5][3] == 'LYS' and Yvalue10[cnt5][4] == sys.argv[4] and Yvalue10[cnt5][2] == 'NZ':
                dist_378 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue10[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue10[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue10[cnt5][8])) ** 2)
                a_378 = [float(Xvalue15[cnt4][5]), float(Yvalue10[cnt5][5]), dist_378]
                # print a
                dist378.append(a_378)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue11)):
            if Yvalue11[cnt5][3] == 'SER' and Yvalue11[cnt5][4] == sys.argv[4] and Yvalue11[cnt5][2] == 'OG':
                dist_379 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue11[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue11[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue11[cnt5][8])) ** 2)
                a_379 = [float(Xvalue15[cnt4][5]), float(Yvalue11[cnt5][5]), dist_379]
                # print a
                dist379.append(a_379)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue12)):
            if Yvalue12[cnt5][3] == 'THR' and Yvalue12[cnt5][4] == sys.argv[4] and Yvalue12[cnt5][2] == 'OG1':
                dist_380 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue12[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue12[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue12[cnt5][8])) ** 2)
                a_380 = [float(Xvalue15[cnt4][5]), float(Yvalue12[cnt5][5]), dist_380]
                # print a
                dist380.append(a_380)



#-backbone to chemical ligand (backbone to side chain) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue13)):
            if Yvalue13[cnt5][3] == 'TRP' and Yvalue13[cnt5][4] == sys.argv[4] and Yvalue13[cnt5][2] == 'NE1':
                dist_381 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue13[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue13[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue13[cnt5][8])) ** 2)
                a_381 = [float(Xvalue15[cnt4][5]), float(Yvalue13[cnt5][5]), dist_381]
                # print a
                dist381.append(a_381)

#-backbone to chemical ligand (backbone to side chain ) x
for cnt4 in range(len(Xvalue15)):
    if Xvalue15[cnt4][3] == 'TYR' and Xvalue15[cnt4][4] == sys.argv[3] and Xvalue15[cnt4][2] == 'O':
        for cnt5 in range(len(Yvalue14)):
            if Yvalue14[cnt5][3] == 'TYR' and Yvalue14[cnt5][4] == sys.argv[4] and Yvalue14[cnt5][2] == 'OH':
                dist_382 = np.sqrt((float(Xvalue15[cnt4][6]) - float(Yvalue14[cnt5][6])) ** 2 + (
                        float(Xvalue15[cnt4][7]) - float(Yvalue14[cnt5][7])) ** 2 + (
                                         float(Xvalue15[cnt4][8]) - float(Yvalue14[cnt5][8])) ** 2)
                a_382 = [float(Xvalue15[cnt4][5]), float(Yvalue14[cnt5][5]), dist_382]
                # print a
                dist382.append(a_382)
#######################################################################################


#output files
#ASNOD1-DONOR
if sys.argv[2] == 'ASN':
    np.savetxt('ASNOD1_ARGNE.csv',dist1,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_ARGNH1.csv',dist2,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_ARGNH2.csv',dist3,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_ASNND2.csv',dist4,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_GLNNE2.csv',dist5,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_HISND1.csv',dist6,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_HISNE2.csv',dist7,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_LYSNZ.csv',dist8,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_SEROG.csv',dist9,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_THROG1.csv',dist10,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_TRPNE1.csv',dist11,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_TYROH.csv',dist12,'%5.2f',delimiter=',')
    np.savetxt('ASNOD1_N12_substrate.csv',dist13,'%5.2f',delimiter=',')
    np.savetxt('ASN_O_mainchain_N_substrate.csv',dist14,'%5.2f',delimiter=',')
    np.savetxt('ASN_O_GLN_N_mainchain.csv',dist15, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_GLU_N_mainchain.csv', dist16, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ASP_N_mainchain.csv', dist17, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_VAL_N_mainchain.csv', dist18, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ILE_N_mainchain.csv', dist19, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_LEU_N_mainchain.csv', dist20, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_MET_N_mainchain.csv', dist21, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_SER_N_mainchain.csv', dist22, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_CYS_N_mainchain.csv', dist23, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_TRP_N_mainchain.csv', dist24, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_TYR_N_mainchain.csv', dist25, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_PRO_N_mainchain.csv', dist26, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ALA_N_mainchain.csv', dist27, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_HIS_N_mainchain.csv', dist28, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ARG_N_mainchain.csv', dist29, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_LYS_N_mainchain.csv', dist30, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_PHE_N_mainchain.csv', dist31, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_THR_N_mainchain.csv', dist32, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ARG_NE_sidechain.csv', dist33, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ARG_NH1_sidechain.csv', dist34, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ARG_NH2_sidechain.csv', dist35, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_ASN_ND2_sidechain.csv', dist36, '%5.2f',delimiter=',' )
    np.savetxt('ASN_O_GLN_NE2_sidechain.csv', dist37, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_HIS_ND1_sidechain.csv', dist38, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_HIS_NE2_sidechain.csv', dist39, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_LYS_NZ_sidechain.csv', dist40, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_SER_OG_sidechain.csv', dist41, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_THR_OG1_sidechain.csv', dist42, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_TRP_NE1_sidechain.csv', dist43, '%5.2f', delimiter=',')
    np.savetxt('ASN_O_TYR_OH_sidechain.csv', dist44, '%5.2f', delimiter=',')


#GLNOE1-DONOR
if sys.argv[2] == 'GLN':
    np.savetxt('GLNOE1_ARGNE.csv',dist45,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_ARGNH1.csv',dist46,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_ARGNH2.csv',dist47,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_ASNND2.csv',dist48,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_GLNNE2.csv',dist49,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_HISND1.csv',dist50,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_HISNE2.csv',dist51,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_LYSNZ.csv',dist52,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_SEROG.csv',dist53,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_THROG1.csv',dist54,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_TRPNE1.csv',dist55,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_TYROH.csv',dist56,'%5.2f',delimiter=',')
    np.savetxt('GLNOE1_N12_substrate.csv',dist57,'%5.2f',delimiter=',')
    np.savetxt('GLN_O_mainchain_N_substrate.csv',dist58,'%5.2f',delimiter=',')
    np.savetxt('GLN_O_GLN_N_mainchain.csv',dist59, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_GLU_N_mainchain.csv', dist60, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ASP_N_mainchain.csv', dist61, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_VAL_N_mainchain.csv', dist62, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ILE_N_mainchain.csv', dist63, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_LEU_N_mainchain.csv', dist64, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_MET_N_mainchain.csv', dist65, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_SER_N_mainchain.csv', dist66, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_CYS_N_mainchain.csv', dist67, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_TRP_N_mainchain.csv', dist68, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_TYR_N_mainchain.csv', dist69, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_PRO_N_mainchain.csv', dist70, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ALA_N_mainchain.csv', dist71, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_HIS_N_mainchain.csv', dist72, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ARG_N_mainchain.csv', dist73, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_LYS_N_mainchain.csv', dist74, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_PHE_N_mainchain.csv', dist75, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_THR_N_mainchain.csv', dist76, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ARG_NE_sidechain.csv', dist77, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ARG_NH1_sidechain.csv', dist78, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ARG_NH2_sidechain.csv', dist79, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_ASN_ND2_sidechain.csv', dist80, '%5.2f',delimiter=',' )
    np.savetxt('GLN_O_GLN_NE2_sidechain.csv', dist81, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_HIS_ND1_sidechain.csv', dist82, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_HIS_NE2_sidechain.csv', dist83, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_LYS_NZ_sidechain.csv', dist84, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_SER_OG_sidechain.csv', dist85, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_THR_OG1_sidechain.csv', dist86, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_TRP_NE1_sidechain.csv', dist87, '%5.2f', delimiter=',')
    np.savetxt('GLN_O_TYR_OH_sidechain.csv', dist88, '%5.2f', delimiter=',')


#GLUOE1-DONOR
#GLUOE2-DONOR
if sys.argv[2] == 'GLU':
    np.savetxt('GLUOE1_ARGNE.csv',dist89,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_ARGNE.csv',dist90,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_ARGNH1.csv',dist91,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_ARGNH1.csv',dist92,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_ARGNH2.csv',dist93,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_ARGNH2.csv',dist94,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_ASNND2.csv',dist95,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_ASNND2.csv',dist96,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_GLNNE2.csv',dist97,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_GLNNE2.csv',dist98,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_HISND1.csv',dist99,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_HISNE2.csv',dist100,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_HISND1.csv',dist101,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_HISNE2.csv',dist102,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_LYSNZ.csv',dist103,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_LYSNZ.csv',dist104,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_SEROG.csv',dist105,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_SEROG.csv',dist106,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_THROG1.csv',dist107,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_THROG1.csv',dist108,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_TRPNE1.csv',dist109,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_TRPNE1.csv',dist110,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_TYROH.csv',dist111,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_TYROH.csv',dist112,'%5.2f',delimiter=',')
    np.savetxt('GLUOE1_N12_substrate.csv',dist113,'%5.2f',delimiter=',')
    np.savetxt('GLUOE2_N12_substrate.csv',dist114,'%5.2f',delimiter=',')
    np.savetxt('GLU_O_mainchain_N_substrate.csv',dist115,'%5.2f',delimiter=',')
    np.savetxt('GLU_O_GLN_N_mainchain.csv',dist116, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_GLU_N_mainchain.csv', dist117, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ASP_N_mainchain.csv', dist118, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_VAL_N_mainchain.csv', dist119, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ILE_N_mainchain.csv', dist120, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_LEU_N_mainchain.csv', dist121, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_MET_N_mainchain.csv', dist122, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_SER_N_mainchain.csv', dist123, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_CYS_N_mainchain.csv', dist124, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_TRP_N_mainchain.csv', dist125, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_TYR_N_mainchain.csv', dist126, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_PRO_N_mainchain.csv', dist127, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ALA_N_mainchain.csv', dist128, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_HIS_N_mainchain.csv', dist129, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ARG_N_mainchain.csv', dist130, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_LYS_N_mainchain.csv', dist131, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_PHE_N_mainchain.csv', dist132, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_THR_N_mainchain.csv', dist133, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ARG_NE_sidechain.csv', dist134, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ARG_NH1_sidechain.csv', dist135, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ARG_NH2_sidechain.csv', dist136, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_ASN_ND2_sidechain.csv', dist137, '%5.2f',delimiter=',' )
    np.savetxt('GLU_O_GLN_NE2_sidechain.csv', dist138, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_HIS_ND1_sidechain.csv', dist139, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_HIS_NE2_sidechain.csv', dist140, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_LYS_NZ_sidechain.csv', dist141, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_SER_OG_sidechain.csv', dist142, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_THR_OG1_sidechain.csv', dist143, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_TRP_NE1_sidechain.csv', dist144, '%5.2f', delimiter=',')
    np.savetxt('GLU_O_TYR_OH_sidechain.csv', dist145, '%5.2f', delimiter=',')

#ASPOD1-DONOR
#ASPOD2-DONOR
if sys.argv[2] == 'ASP':
    np.savetxt('ASPOD1_ARGNE.csv',dist146,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_ARGNE.csv',dist147,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_ARGNH1.csv',dist148,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_ARGNH1.csv',dist149,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_ARGNH2.csv',dist150,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_ARGNH2.csv',dist151,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_ASNND2.csv',dist152,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_ASNND2.csv',dist153,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_GLNNE2.csv',dist154,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_GLNNE2.csv',dist155,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_HISND1.csv',dist156,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_HISNE2.csv',dist157,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_HISND1.csv',dist158,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_HISNE2.csv',dist159,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_LYSNZ.csv',dist160,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_LYSNZ.csv',dist161,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_SEROG.csv',dist162,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_SEROG.csv',dist163,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_THROG1.csv',dist164,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_THROG1.csv',dist165,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_TRPNE1.csv',dist166,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_TRPNE1.csv',dist167,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_TYROH.csv',dist168,'%5.2f',delimiter=',')
    np.savetxt('ASPOD2_TYROH.csv',dist169,'%5.2f',delimiter=',')
    np.savetxt('ASPOD1_N12_substrate.csv',dist170,'%5.2f',delimiter=',')
    np.savetxt('ASP_O_mainchain_N_substrate.csv',dist171,'%5.2f',delimiter=',')
    np.savetxt('ASP_O_GLN_N_mainchain.csv',dist172, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_GLU_N_mainchain.csv', dist173, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ASP_N_mainchain.csv', dist174, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_VAL_N_mainchain.csv', dist175, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ILE_N_mainchain.csv', dist176, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_LEU_N_mainchain.csv', dist177, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_MET_N_mainchain.csv', dist178, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_SER_N_mainchain.csv', dist179, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_CYS_N_mainchain.csv', dist180, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_TRP_N_mainchain.csv', dist181, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_TYR_N_mainchain.csv', dist182, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_PRO_N_mainchain.csv', dist183, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ALA_N_mainchain.csv', dist184, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_HIS_N_mainchain.csv', dist185, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ARG_N_mainchain.csv', dist186, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_LYS_N_mainchain.csv', dist187, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_PHE_N_mainchain.csv', dist188, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_THR_N_mainchain.csv', dist189, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ARG_NE_sidechain.csv', dist190, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ARG_NH1_sidechain.csv', dist191, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ARG_NH2_sidechain.csv', dist192, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_ASN_ND2_sidechain.csv', dist193, '%5.2f',delimiter=',' )
    np.savetxt('ASP_O_GLN_NE2_sidechain.csv', dist194, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_HIS_ND1_sidechain.csv', dist195, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_HIS_NE2_sidechain.csv', dist196, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_LYS_NZ_sidechain.csv', dist197, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_SER_OG_sidechain.csv', dist198, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_THR_OG1_sidechain.csv', dist199, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_TRP_NE1_sidechain.csv', dist200, '%5.2f', delimiter=',')
    np.savetxt('ASP_O_TYR_OH_sidechain.csv', dist201, '%5.2f', delimiter=',')


#HISNE2-DONOR
if sys.argv[2] == 'HSD':
    np.savetxt('HISND1_ARGNE.csv',dist202,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_ARGNE.csv',dist203,'%5.2f',delimiter=',')
    np.savetxt('HISND1_ARGNH1.csv',dist204,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_ARGNH1.csv',dist205,'%5.2f',delimiter=',')
    np.savetxt('HISND1_ARGNH2.csv',dist206,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_ARGNH2.csv',dist207,'%5.2f',delimiter=',')
    np.savetxt('HISND1_ASNND2.csv',dist208,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_ASNND2.csv',dist209,'%5.2f',delimiter=',')
    np.savetxt('HISND1_GLNNE2.csv',dist210,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_GLNNE2.csv',dist211,'%5.2f',delimiter=',')
    np.savetxt('HISND1_HISND1.csv',dist212,'%5.2f',delimiter=',')
    np.savetxt('HISND1_HISNE2.csv',dist213,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_HISND1.csv',dist214,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_HISNE2.csv',dist215,'%5.2f',delimiter=',')
    np.savetxt('HISND1_LYSNZ.csv',dist216,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_LYSNZ.csv',dist217,'%5.2f',delimiter=',')
    np.savetxt('HISND1_SEROG.csv',dist218,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_SEROG.csv',dist219,'%5.2f',delimiter=',')
    np.savetxt('HISND1_THROG1.csv',dist220,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_THROG1.csv',dist221,'%5.2f',delimiter=',')
    np.savetxt('HISND1_TRPNE1.csv',dist222,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_TRPNE1.csv',dist223,'%5.2f',delimiter=',')
    np.savetxt('HISND1_TYROH.csv',dist224,'%5.2f',delimiter=',')
    np.savetxt('HISNE2_TYROH.csv',dist225,'%5.2f',delimiter=',')
    np.savetxt('HISND1_N12_substrate.csv',dist226,'%5.2f',delimiter=',')
    np.savetxt('HIS_O_N12_substrate.csv',dist227,'%5.2f',delimiter=',')
    np.savetxt('HIS_O_HIS_N_mainchain', dist228,'%5.2f', delimiter=',')
    np.savetxt('HIS_O_VAL_N_mainchain', dist229, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ILE_N_mainchain', dist230, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_LEU_N_mainchain', dist231, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_MET_N_mainchain', dist232, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_SER_N_mainchain', dist233, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_CYS_N_mainchain', dist234, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_TRP_N_mainchain', dist235, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_TYR_N_mainchain', dist236, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_PRO_N_mainchain', dist237, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ALA_N_mainchain', dist238, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_HIS_N_mainchain', dist239, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ARG_N_mainchain', dist240, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_LYS_N_mainchain', dist241, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_PHE_N_mainchain', dist242, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_THR_N_mainchain', dist243, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ARG_NE_sidechain', dist244, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ARG_NH1_sidechain', dist245, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ARG_NH2_sidechain', dist246, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_ASN_ND2_sidechain', dist247, '%5.2f',delimiter=',' )
    np.savetxt('HIS_O_GLN_NE2_sidechain', dist248, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_HIS_ND1_sidechain', dist249, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_HIS_NE2_sidechain', dist250, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_LYS_NZ_sidechain', dist251, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_SER_OG_sidechain', dist252, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_THR_OG1_sidechain', dist253, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_TRP_NE1_sidechain', dist254, '%5.2f', delimiter=',')
    np.savetxt('HIS_O_TYR_OH_sidechain', dist255, '%5.2f', delimiter=',')


# ACCEPTOR-DONOR SER
if sys.argv[2] == 'SER':
    np.savetxt('SEROG_ARGNE.csv',dist256,'%5.2f',delimiter=',')
    np.savetxt('SEROG_ARGNH1.csv',dis257,'%5.2f',delimiter=',')
    np.savetxt('SEROG_ARGNH2.csv',dis258,'%5.2f',delimiter=',')
    np.savetxt('SEROG_ASNND2.csv',dist259,'%5.2f',delimiter=',')
    np.savetxt('SEROG_GLNNE2.csv',dist260,'%5.2f',delimiter=',')
    np.savetxt('SEROG_HISND1.csv',dist261,'%5.2f',delimiter=',')
    np.savetxt('SEROG_HISNE2.csv',dist262,'%5.2f',delimiter=',')
    np.savetxt('SEROG_LYSNZ.csv',dist263,'%5.2f',delimiter=',')
    np.savetxt('SEROG_SEROG.csv',dist264,'%5.2f',delimiter=',')
    np.savetxt('SEROG_THROG1.csv',dist265,'%5.2f',delimiter=',')
    np.savetxt('SEROG_TRPNE1.csv',dist266,'%5.2f',delimiter=',')
    np.savetxt('SEROG_TYROH.csv',dist267,'%5.2f',delimiter=',')
    np.savetxt('SEROE1_N12_substrate.csv',dist268,'%5.2f',delimiter=',')
    np.savetxt('SER_O_mainchain_N_substrate.csv',dist269,'%5.2f',delimiter=',')
    np.savetxt('SER_O_VAL_N_mainchain.csv', dist270, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ILE_N_mainchain.csv', dist271, '%5.2f', delimiter=',')
    np.savetxt('SER_O_LEU_N_mainchain.csv', dist272, '%5.2f', delimiter=',')
    np.savetxt('SER_O_MET_N_mainchain.csv', dist273, '%5.2f', delimiter=',')
    np.savetxt('SER_O_SER_N_mainchain.csv', dist274, '%5.2f', delimiter=',')
    np.savetxt('SER_O_CYS_N_mainchain.csv', dist275, '%5.2f', delimiter=',')
    np.savetxt('SER_O_TRP_N_mainchain.csv', dist276, '%5.2f', delimiter=',')
    np.savetxt('SER_O_TYR_N_mainchain.csv', dist277, '%5.2f', delimiter=',')
    np.savetxt('SER_O_PRO_N_mainchain.csv', dist278, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ALA_N_mainchain.csv', dist279, '%5.2f', delimiter=',')
    np.savetxt('SER_O_HIS_N_mainchain.csv', dist280, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ARG_N_mainchain.csv', dist281, '%5.2f', delimiter=',')
    np.savetxt('SER_O_LYS_N_mainchain.csv', dist282, '%5.2f', delimiter=',')
    np.savetxt('SER_O_PHE_N_mainchain.csv', dist283, '%5.2f', delimiter=',')
    np.savetxt('SER_O_THR_N_mainchain.csv', dist284, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ARG_NE_sidechain.csv', dist285, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ARG_NH1_sidechain.csv', dist286, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ARG_NH2_sidechain.csv', dist287, '%5.2f', delimiter=',')
    np.savetxt('SER_O_ASN_ND2_sidechain.csv', dist288, '%5.2f',delimiter=',' )
    np.savetxt('SER_O_GLN_NE2_sidechain.csv', dist289, '%5.2f', delimiter=',')
    np.savetxt('SER_O_HIS_ND1_sidechain.csv', dist290, '%5.2f', delimiter=',')
    np.savetxt('SER_O_HIS_NE2_sidechain.csv', dist291, '%5.2f', delimiter=',')
    np.savetxt('SER_O_LYS_NZ_sidechain.csv', dist292, '%5.2f', delimiter=',')
    np.savetxt('SER_O_SER_OG_sidechain.csv', dist293, '%5.2f', delimiter=',')
    np.savetxt('SER_O_THR_OG1_sidechain.csv', dist294, '%5.2f', delimiter=',')
    np.savetxt('SER_O_TRP_NE1_sidechain.csv', dist295, '%5.2f', delimiter=',')
    np.savetxt('SER_O_TYR_OH_sidechain.csv', dist296, '%5.2f', delimiter=',')


# ACCEPTOR-DONOR THR
if sys.argv[2] == 'THR':
    np.savetxt('THROG1_ARGNE.csv',dist297,'%5.2f',delimiter=',')
    np.savetxt('THROG1_ARGNH1.csv',dist298,'%5.2f',delimiter=',')
    np.savetxt('THROG1_ARGNH2.csv',dist299,'%5.2f',delimiter=',')
    np.savetxt('THROG1_ASNND2.csv',dist300,'%5.2f',delimiter=',')
    np.savetxt('THROG1_GLNNE2.csv',dist301,'%5.2f',delimiter=',')
    np.savetxt('THROG1_HISND1.csv',dist302,'%5.2f',delimiter=',')
    np.savetxt('THROG1_HISNE2.csv',dist303,'%5.2f',delimiter=',')
    np.savetxt('THROG1_LYSNZ.csv',dist304,'%5.2f',delimiter=',')
    np.savetxt('THROG1_SEROG.csv',dist305,'%5.2f',delimiter=',')
    np.savetxt('THROG1_THROG1.csv',dist306,'%5.2f',delimiter=',')
    np.savetxt('THROG1_TRPNE1.csv',dist307,'%5.2f',delimiter=',')
    np.savetxt('THROG1_TYROH.csv',dist308,'%5.2f',delimiter=',')
    np.savetxt('THROE1_N12_substrate.csv',dist309,'%5.2f',delimiter=',')
    np.savetxt('THR_O_mainchain_N_substrate.csv',dist310,'%5.2f',delimiter=',')
    np.savetxt('THR_O_GLN_N_mainchain.csv',dist311, '%5.2f', delimiter=',')
    np.savetxt('THR_O_GLU_N_mainchain.csv', dist312, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ASP_N_mainchain.csv', dist313, '%5.2f', delimiter=',')
    np.savetxt('THR_O_VAL_N_mainchain.csv', dist314, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ILE_N_mainchain.csv', dist315, '%5.2f', delimiter=',')
    np.savetxt('THR_O_LEU_N_mainchain.csv', dist316, '%5.2f', delimiter=',')
    np.savetxt('THR_O_MET_N_mainchain.csv', dist317, '%5.2f', delimiter=',')
    np.savetxt('THR_O_SER_N_mainchain.csv', dist318, '%5.2f', delimiter=',')
    np.savetxt('THR_O_CYS_N_mainchain.csv', dist319, '%5.2f', delimiter=',')
    np.savetxt('THR_O_TRP_N_mainchain.csv', dist320, '%5.2f', delimiter=',')
    np.savetxt('THR_O_TYR_N_mainchain.csv', dist321, '%5.2f', delimiter=',')
    np.savetxt('THR_O_PRO_N_mainchain.csv', dist322, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ALA_N_mainchain.csv', dist323, '%5.2f', delimiter=',')
    np.savetxt('THR_O_HIS_N_mainchain.csv', dist324, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ARG_N_mainchain.csv', dist325, '%5.2f', delimiter=',')
    np.savetxt('THR_O_LYS_N_mainchain.csv', dist326, '%5.2f', delimiter=',')
    np.savetxt('THR_O_PHE_N_mainchain.csv', dist327, '%5.2f', delimiter=',')
    np.savetxt('THR_O_THR_N_mainchain.csv', dist328, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ARG_NE_sidechain.csv', dist329, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ARG_NH1_sidechain.csv', dist330, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ARG_NH2_sidechain.csv', dist331, '%5.2f', delimiter=',')
    np.savetxt('THR_O_ASN_ND2_sidechain.csv', dist332, '%5.2f',delimiter=',' )
    np.savetxt('THR_O_GLN_NE2_sidechain.csv', dist332, '%5.2f', delimiter=',')
    np.savetxt('THR_O_HIS_ND1_sidechain.csv', dist333, '%5.2f', delimiter=',')
    np.savetxt('THR_O_HIS_NE2_sidechain.csv', dist334, '%5.2f', delimiter=',')
    np.savetxt('THR_O_LYS_NZ_sidechain.csv', dist335, '%5.2f', delimiter=',')
    np.savetxt('THR_O_SER_OG_sidechain.csv', dist336, '%5.2f', delimiter=',')
    np.savetxt('THR_O_THR_OG1_sidechain.csv', dist337, '%5.2f', delimiter=',')
    np.savetxt('THR_O_TRP_NE1_sidechain.csv', dist338, '%5.2f', delimiter=',')
    np.savetxt('THR_O_TYR_OH_sidechain.csv', dist339, '%5.2f', delimiter=',')


# ACCEPTOR_DONOR TYR
if sys.argv[2] == 'TYR':
    np.savetxt('TYROH_ARGNE.csv',dist340,'%5.2f',delimiter=',')
    np.savetxt('TYROH_ARGNH1.csv',dist341,'%5.2f',delimiter=',')
    np.savetxt('TYROH_ARGNH2.csv',dist342,'%5.2f',delimiter=',')
    np.savetxt('TYROH_ASNND2.csv',dist343,'%5.2f',delimiter=',')
    np.savetxt('TYROH_GLNNE2.csv',dist344,'%5.2f',delimiter=',')
    np.savetxt('TYROH_HISND1.csv',dist345,'%5.2f',delimiter=',')
    np.savetxt('TYROH_HISNE2.csv',dist346,'%5.2f',delimiter=',')
    np.savetxt('TYROH_LYSNZ.csv',dist347,'%5.2f',delimiter=',')
    np.savetxt('TYROH_SEROH.csv',dist348,'%5.2f',delimiter=',')
    np.savetxt('TYROH_THROG1.csv',dist349,'%5.2f',delimiter=',')
    np.savetxt('TYROH_TRPNE1.csv',dist350,'%5.2f',delimiter=',')
    np.savetxt('TYROH_TYROH.csv',dist351,'%5.2f',delimiter=',')
    np.savetxt('TYROE1_N12_substrate.csv',dist352,'%5.2f',delimiter=',')
    np.savetxt('TYR_O_mainchain_N_substrate.csv',dist352,'%5.2f',delimiter=',')
    np.savetxt('TYR_O_GLN_N_mainchain.csv',dist353, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_GLU_N_mainchain.csv', dist354, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ASP_N_mainchain.csv', dist355, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_VAL_N_mainchain.csv', dist356, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ILE_N_mainchain.csv', dist357, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_LEU_N_mainchain.csv', dist358, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_MET_N_mainchain.csv', dist359, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_SER_N_mainchain.csv', dist360, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_CYS_N_mainchain.csv', dist361, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_TRP_N_mainchain.csv', dist362, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_TYR_N_mainchain.csv', dist363, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_PRO_N_mainchain.csv', dist364, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ALA_N_mainchain.csv', dist365, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_HIS_N_mainchain.csv', dist366, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ARG_N_mainchain.csv', dist367, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_LYS_N_mainchain.csv', dist368, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_PHE_N_mainchain.csv', dist369, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_THR_N_mainchain.csv', dist370, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ARG_NE_sidechain.csv', dist371, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ARG_NH1_sidechain.csv', dist372, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ARG_NH2_sidechain.csv', dist373, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_ASN_ND2_sidechain.csv', dist374, '%5.2f',delimiter=',' )
    np.savetxt('TYR_O_GLN_NE2_sidechain.csv', dist375, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_HIS_ND1_sidechain.csv', dist376, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_HIS_NE2_sidechain.csv', dist377, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_LYS_NZ_sidechain.csv', dist378, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_SER_OG_sidechain.csv', dist379, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_THR_OG1_sidechain.csv', dist380, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_TRP_NE1_sidechain.csv', dist381, '%5.2f', delimiter=',')
    np.savetxt('TYR_O_TYR_OH_sidechain.csv', dist382, '%5.2f', delimiter=',')
