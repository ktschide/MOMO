function [random_vector] = random_LSareas(myfun,powermin,powermax,areas)


x = logspace(powermin,powermax,100);

mypdf = myfun(x);

random_vector = randsample(x,areas,true,mypdf);


