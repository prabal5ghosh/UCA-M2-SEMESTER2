b = 1;
a = [1, 0.5]
n = 0:19;
x = cos(n*pi/2)
y = filter(b,a,x)

figure

stem(n,x)
hold on
stem(n,y)


% theoritical expression
yt = 2/sqrt(5)*cos(n*pi/2 - atan(0.5))

stem(n,yt)

grid(gce, 'FontSize', fontsize )

xlabel('sample number, n')
ylabel('y(n)')
