function output = circtukey3D(dimz,dimy,dimx,lev,row,col,filterwidth)

% 3D Tukey filter

domain = 512;

base = zeros(domain,domain,domain);

tukey1 = tukeywin(domain,filterwidth);
tukey1 = tukey1(domain/2+1:domain);

shiftz = (lev-dimz/2)*domain/dimz;
shifty = (row-dimy/2)*domain/dimy;
shiftx = (col-dimx/2)*domain/dimx;

z = linspace(-domain/2, domain/2, domain);
y = linspace(-domain/2, domain/2, domain);
x = linspace(-domain/2, domain/2, domain);

for i=1:domain

    for j=1:domain
    
        for k = 1:domain
        
            rad = round(sqrt((shiftx-x(i))^2 + (shifty-y(j))^2 + (shiftz-z(k))^2));
            
            if (rad <= domain/2) && (rad > 0)
                
              base(k,j,i) = tukey1(rad);
              
            end
            
        end
        
    end
    
end

output = imresize3(base,[dimz dimy dimx]);


end