function output = circtukey2D(dimy,dimx,filterwidth)

% 2D Tukey filter

domain = 512;
base = zeros(domain,domain);

tukey1 = tukeywin(domain,filterwidth);
tukey1 = tukey1(domain/2+1:domain);

x = linspace(-domain/2, domain/2, domain);
y = linspace(-domain/2, domain/2, domain);

for i=1:domain
    
    for j=1:domain
        
        if (round(sqrt(x(i)^2 + y(j)^2)) <= domain/2)
            
            base(i,j) = tukey1(round(sqrt(x(i)^2 + y(j)^2)));
            
        end
        
    end
    
end

output = imresize(base,[dimy dimx]);

end