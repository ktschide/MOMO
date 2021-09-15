function [Locyarray, Locxarray, la] = fitlandslide(la, cellsize_m2, scape, xloc,yloc);

            [row, col] = size(scape);
            disturbed_cells=la/(cellsize_m2);
            sqcells=round(sqrt(disturbed_cells));
            
            %%fit size to array and location
            ydir=sqcells;
            xdir=sqcells;
             
            if (ydir-1+yloc)>row
                ydir=row-yloc;
            end
            
            if (xdir-1+xloc)>col
                xdir=col-xloc;
            end
           
          
            yprop = ydir-1+yloc;
            Lsrows = scape(yloc:yprop,xloc);   
            
            xprop = xdir-1+xloc;
            Lscols = scape(yloc, xloc:xprop);   
            
            Locyarray=(yloc:1:(yloc+ydir-1));
            Locxarray=(xloc:1:(xloc+xdir-1));
            
            
            if isnan(sum(Lsrows))
                [yedge,~] = find(isnan(Lsrows));
                minnany = 1:min(yedge)-1;
                Locyarray = (yloc+(minnany-1))';
            end
            if isnan(sum(Lscols))
                [~,xedge] = find(isnan(Lscols));
                minnanx = 1:min(xedge)-1;
                Locxarray = (xloc+(minnanx-1))';
            end
            [nanrow, nancol]  = find(isnan(scape(Locyarray, Locxarray)));
            edge_effect = length(nanrow);
           
            la = ((length(Locyarray))*(length(Locxarray))-edge_effect)*cellsize_m2;
            
            
end