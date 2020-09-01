function z = Sphere(Irow,maxRowC)
tempFunc=0;
for i=1:9
    tempFunc=tempFunc+(maxRowC-Irow(i));
end
z=tempFunc;
%z=sum(panelMax)+(10/((maxRowC-Irow1)+(maxRowC-Irow2)+(maxRowC-Irow3)));

end

%+(10./(maxRowC-Irow))+(10./(maxRowC-Irow))+(10.*panelMax)