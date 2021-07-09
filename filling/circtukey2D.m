function output = circtukey2D(dimy,dimx,row,col,filterwidth)

% 2D Tukey filter

domain = 512;
base = zeros(domain,domain);

tukey1 = tukeywin(domain,filterwidth);
tukey1 = tukey1(domain/2+1:domain);

shifty = (row-dimy/2)*domain/dimy;
shiftx = (col-dimx/2)*domain/dimx;

y = linspace(-domain/2, domain/2, domain);
x = linspace(-domain/2, domain/2, domain);

for i=1:domain
    
    for j=1:domain
        
        rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2)); 
        
        if (rad <= domain/2) && (rad > 0)
            
            base(j,i) = tukey1(rad);
            
        end
        
    end
    
end

output = imresize(base,[dimy dimx]);

end
