% This is the file parabola.m
% Calculate squares of array
% Input: x      Return: y
function[y] = parabola(x)
    [rows,cols] = size(x);
    for i = 1:cols
        y(i) = x(i) ^2;
    end
end