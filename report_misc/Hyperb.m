syms x y z
c=1;
fimplicit3(sqrt((x-c)^2 + y^2 + z^2)-sqrt((x+c)^2 + y^2 + z^2)-c)
hold on;
fimplicit3(sqrt((x-c)^2 + y^2 + z^2)-sqrt((x+c)^2 + y^2 + z^2)+c)