macro qx(ix, iy) esc(:( (C[$ix+1, $iy] - C[$ix, $iy])  *-D /dx^2 )) end
macro qy(ix, iy) esc(:( (C[$ix, $iy+1] - C[$ix, $iy])  *-D /dy^2 )) end

macro fcx(ix, iy) esc(:( 0.5 /dx *(
    C[$ix+1, $iy] + C[$ix, $iy]) *Uf[$ix, $iy] )) end
macro fcy(ix, iy) esc(:( 0.5 /dy *(
    C[$ix, $iy+1] + C[$ix, $iy]) *Vf[$ix, $iy] )) end

function step_transport!(C2, C, Uf, Vf, DeltaX, DeltaT) # --------------------
    # nG = 1
    # ImaAll, JmaAll = size(C)
	# Ifi = nG+1;     Ila = ImaAll - nG
	# Jfi = nG+1;     Jla = JmaAll - nG

    D = 1e-5
    dx = dy = DeltaX
    dt = DeltaT

    for j in Jfi:Jla
        for i in Ifi:Ila
            @inbounds C2[i, j] = C[i, j] + dt * ( 
                - @fcx(i, j) + @fcx(i-1, j) - @fcy(i, j) + @fcx(i, j-1)
                - @qx(i, j) + @qx(i-1, j) - @qy(i, j) + @qy(i, j-1)
             )
        end
    end
    return nothing
end # ----------------------------------------------------------------

