function output = circtukey3D(dimz,dimy,dimx,filterwidth)

% 3D Tukey filter

domain = 512;

base = zeros(domain,domain,domain);

tukey1 = tukeywin(domain,filterwidth);
tukey1 = tukey1(domain/2+1:domain);

x = linspace(-domain/2, domain/2, domain);
y = linspace(-domain/2, domain/2, domain);
z = linspace(-domain/2, domain/2, domain);

for i=1:domain

    for j=1:domain
    
        for k = 1:domain
        
            if (round(sqrt(x(i)^2 + y(j)^2 + z(k)^2)) <= domain/2)
                
              base(i,j,k) = tukey1(round(sqrt(x(i)^2 + y(j)^2 + z(k)^2)));
              
            end
            
        end
        
    end
    
end

output = imresize3(base,[dimz dimy dimx]);

end