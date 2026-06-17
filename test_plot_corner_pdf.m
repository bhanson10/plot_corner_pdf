figure; hold on; 
d = 5; mu = randn(1,d); P = randn(d,d); P = P' * P; N = 20; x = zeros(d, N); 
for i = 1:d, x(i,:) = linspace(-3 * sqrt(P(i,i)) + mu(i),3 * sqrt(P(i,i)) + mu(i),N); end
[X1,X2,X3,X4,X5] = ndgrid(x(1,:)',x(2,:)',x(3,:)',x(4,:)',x(5,:)'); 
X = [X1(:) X2(:) X3(:) X4(:) X5(:)]; 
p.color = 'r'; p.type = "grid"; 
P = mvnpdf(X, mu, P); plot_corner_pdf(X, 'P', P, 'p', p);