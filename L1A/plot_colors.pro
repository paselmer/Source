FUNCTION plot_colors;, cindex_array

;INPUT:  None. Simply set a variable equal to this function with a
;        blank call
;OUTPUT: 1x10 longword integer array with color indices [cindex_array]

;DESCRIPTION:
;This function returns an array filled with color index values.
;These indices can be fed to the "color" keyword in the
;plot routine.
;The function fills the array with color indices calculated from
;the hard-coded RGB colors in the "colors" array below.
;
;Written by Patrick Selmer
;Modified - 5/12/13

ncolors=161
cindex_array=lindgen(ncolors)

;different colors
colors=lindgen(3,ncolors)
colors[0:2,0]=[255,105,180] ;hot pink
colors[0:2,1]=[0,255,0]     ;lime
colors[0:2,2]=[32,178,170]  ;light sea green
colors[0:2,3]=[0,0,255]     ;blue
colors[0:2,4]=[210,105,30]  ;chocolate
colors[0:2,5]=[138,43,226]  ;blue violet
colors[0:2,6]=[255,255,0]   ;yellow
colors[0:2,7]=[222,184,135] ;burlywood
colors[0:2,8]=[128,0,0]     ;maroon
colors[0:2,9]=[0,191,255]   ;deep sky blue
colors[0:2,10]=[198,99,231] ;lavender
colors[0:2,11]=[112,128,144];slate gray
colors[0:2,12]=[255,215,0]  ;gold
colors[0:2,13]=[205,92,92]  ;indian red
colors[0:2,14]=[255,0,0]    ;red
colors[0:2,15]=[255,140,0]  ;dark orange
colors[0:2,16]=[238,203,173];peach puff 2
colors[0:2,17]=[238,213,210];misty rose 2
colors[0:2,18]=[142,229,238];cadet blue 2
colors[0:2,19]=[0,139,0]    ;green 4
colors[0:2,20]=[255,236,139];light goldenrod 1
colors[0:2,21]=[139,058,058];indian red 4
colors[0:2,22]=[144,238,144];light green
;Filler colors ---------------
colors[0:2,23]=[255,105,180] ;hot pink
colors[0:2,24]=[0,255,0]     ;lime
colors[0:2,25]=[32,178,170]  ;light sea green
colors[0:2,26]=[0,0,255]     ;blue
colors[0:2,27]=[210,105,30]  ;chocolate
colors[0:2,28]=[138,43,226]  ;blue violet
colors[0:2,29]=[255,255,0]   ;yellow
colors[0:2,30]=[222,184,135] ;burlywood
colors[0:2,31]=[128,0,0]     ;maroon
colors[0:2,32]=[0,191,255]   ;deep sky blue
colors[0:2,33]=[198,99,231] ;lavender
colors[0:2,34]=[112,128,144];slate gray
colors[0:2,35]=[255,215,0]  ;gold
colors[0:2,36]=[205,92,92]  ;indian red
colors[0:2,37]=[255,0,0]    ;red
colors[0:2,38]=[255,140,0]  ;dark orange
colors[0:2,39]=[238,203,173];peach puff 2
colors[0:2,40]=[238,213,210];misty rose 2
colors[0:2,41]=[142,229,238];cadet blue 2
colors[0:2,42]=[0,139,0]    ;green 4
colors[0:2,43]=[255,236,139];light goldenrod 1
colors[0:2,44]=[139,058,058];indian red 4
colors[0:2,45]=[144,238,144];light greencolors[0:2,0]=[255,105,180] ;hot pink
colors[0:2,46]=[0,255,0]     ;lime
colors[0:2,47]=[32,178,170]  ;light sea green
colors[0:2,48]=[0,0,255]     ;blue
colors[0:2,49]=[210,105,30]  ;chocolate
colors[0:2,50]=[138,43,226]  ;blue violet
colors[0:2,51]=[255,255,0]   ;yellow
colors[0:2,52]=[222,184,135] ;burlywood
colors[0:2,53]=[128,0,0]     ;maroon
colors[0:2,54]=[0,191,255]   ;deep sky blue
colors[0:2,55]=[198,99,231] ;lavender
colors[0:2,56]=[112,128,144];slate gray
colors[0:2,57]=[255,215,0]  ;gold
colors[0:2,58]=[205,92,92]  ;indian red
colors[0:2,59]=[255,0,0]    ;red
colors[0:2,60]=[255,140,0]  ;dark orange
colors[0:2,61]=[238,203,173];peach puff 2
colors[0:2,62]=[238,213,210];misty rose 2
colors[0:2,63]=[142,229,238];cadet blue 2
colors[0:2,64]=[0,139,0]    ;green 4
colors[0:2,65]=[255,236,139];light goldenrod 1
colors[0:2,66]=[139,058,058];indian red 4
colors[0:2,67]=[144,238,144];light green
colors[0:2,68]=[255,105,180] ;hot pink
colors[0:2,69]=[0,255,0]     ;lime
colors[0:2,70]=[32,178,170]  ;light sea green
colors[0:2,71]=[0,0,255]     ;blue
colors[0:2,72]=[210,105,30]  ;chocolate
colors[0:2,73]=[138,43,226]  ;blue violet
colors[0:2,74]=[255,255,0]   ;yellow
colors[0:2,75]=[222,184,135] ;burlywood
colors[0:2,76]=[128,0,0]     ;maroon
colors[0:2,77]=[0,191,255]   ;deep sky blue
colors[0:2,78]=[198,99,231] ;lavender
colors[0:2,79]=[112,128,144];slate gray
colors[0:2,80]=[255,215,0]  ;gold
colors[0:2,81]=[205,92,92]  ;indian red
colors[0:2,82]=[255,0,0]    ;red
colors[0:2,83]=[255,140,0]  ;dark orange
colors[0:2,84]=[238,203,173];peach puff 2
colors[0:2,85]=[238,213,210];misty rose 2
colors[0:2,86]=[142,229,238];cadet blue 2
colors[0:2,87]=[0,139,0]    ;green 4
colors[0:2,88]=[255,236,139];light goldenrod 1
colors[0:2,157]=[255,236,139];light goldenrod 1
colors[0:2,89]=[139,058,058];indian red 4
colors[0:2,90]=[144,238,144];light green
colors[0:2,91]=[255,105,180] ;hot pink
colors[0:2,92]=[0,255,0]     ;lime
colors[0:2,93]=[32,178,170]  ;light sea green
colors[0:2,94]=[0,0,255]     ;blue
colors[0:2,95]=[210,105,30]  ;chocolate
colors[0:2,96]=[138,43,226]  ;blue violet
colors[0:2,97]=[255,255,0]   ;yellow
colors[0:2,98]=[222,184,135] ;burlywood
colors[0:2,99]=[128,0,0]     ;maroon
colors[0:2,111]=[255,236,139];light goldenrod 1
colors[0:2,112]=[139,058,058];indian red 4
colors[0:2,113]=[144,238,144];light green
colors[0:2,114]=[255,105,180] ;hot pink
colors[0:2,115]=[0,255,0]     ;lime
colors[0:2,116]=[32,178,170]  ;light sea green
colors[0:2,117]=[0,0,255]     ;blue
colors[0:2,118]=[210,105,30]  ;chocolate
colors[0:2,119]=[138,43,226]  ;blue violet
colors[0:2,120]=[255,255,0]   ;yellow
colors[0:2,121]=[222,184,135] ;burlywood
colors[0:2,122]=[128,0,0]     ;maroon
colors[0:2,123]=[0,191,255]   ;deep sky blue
colors[0:2,124]=[198,99,231] ;lavender
colors[0:2,125]=[112,128,144];slate gray
colors[0:2,100]=[0,191,255]   ;deep sky blue
colors[0:2,101]=[198,99,231] ;lavender
colors[0:2,102]=[112,128,144];slate gray
colors[0:2,103]=[255,215,0]  ;gold
colors[0:2,104]=[205,92,92]  ;indian red
colors[0:2,105]=[255,0,0]    ;red
colors[0:2,106]=[255,140,0]  ;dark orange
colors[0:2,107]=[238,203,173];peach puff 2
colors[0:2,108]=[238,213,210];misty rose 2
colors[0:2,109]=[142,229,238];cadet blue 2
colors[0:2,110]=[0,139,0]    ;green 4
colors[0:2,126]=[255,215,0]  ;gold
colors[0:2,127]=[205,92,92]  ;indian red
colors[0:2,128]=[255,0,0]    ;red
colors[0:2,129]=[255,140,0]  ;dark orange
colors[0:2,130]=[238,203,173];peach puff 2
colors[0:2,152]=[255,140,0]  ;dark orange
colors[0:2,153]=[238,203,173];peach puff 2
colors[0:2,154]=[238,213,210];misty rose 2
colors[0:2,155]=[142,229,238];cadet blue 2
colors[0:2,156]=[0,139,0]    ;green 4
colors[0:2,158]=[139,058,058];indian red 4
colors[0:2,159]=[144,238,144];light green
colors[0:2,160]=[255,255,255] ;white
colors[0:2,131]=[238,213,210];misty rose 2
colors[0:2,132]=[142,229,238];cadet blue 2
colors[0:2,133]=[0,139,0]    ;green 4
colors[0:2,134]=[255,236,139];light goldenrod 1
colors[0:2,135]=[139,058,058];indian red 4
colors[0:2,136]=[144,238,144];light green
colors[0:2,137]=[255,255,255] ;white
colors[0:2,138]=[0,255,0]     ;lime
colors[0:2,139]=[32,178,170]  ;light sea green
colors[0:2,140]=[0,0,255]     ;blue
colors[0:2,141]=[210,105,30]  ;chocolate
colors[0:2,142]=[138,43,226]  ;blue violet
colors[0:2,143]=[255,255,0]   ;yellow
colors[0:2,144]=[222,184,135] ;burlywood
colors[0:2,145]=[128,0,0]     ;maroon
colors[0:2,146]=[0,191,255]   ;deep sky blue
colors[0:2,147]=[198,99,231] ;lavender
colors[0:2,148]=[112,128,144];slate gray
colors[0:2,149]=[255,215,0]  ;gold
colors[0:2,150]=[205,92,92]  ;indian red
colors[0:2,151]=[255,0,0]    ;red
;calculate color index value
for i=0,ncolors-1 do begin
	cindex_array[i] = colors[0,i] + 256L * (colors[1,i] + 256L * colors[2,i])
endfor

return,cindex_array

end
