function [myscaledcolormap]=scaledcolormap(cmaplength,colorcode,minColor,maxColor)
scaledcolor=linspace(minColor,maxColor,cmaplength);
myscaledcolormap=[scaledcolor*colorcode(1);scaledcolor*colorcode(2);scaledcolor*colorcode(3)];
myscaledcolormap=myscaledcolormap';
end
