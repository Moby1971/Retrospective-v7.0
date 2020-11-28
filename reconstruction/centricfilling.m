function scheme = centricfilling(no_views_2)

for g=1:round(no_views_2/2)
    
    ord2(2*g-1)=no_views_2/2+g;
    ord2(2*g)=no_views_2/2-g+1;
    
end

scheme = round(ord2);

end