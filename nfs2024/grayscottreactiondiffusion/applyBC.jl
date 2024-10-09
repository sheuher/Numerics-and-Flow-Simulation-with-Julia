function applyBC!(Phi)
    """
        apply periodic boundary condition
    """
    nG = 1
    ImaAll, JmaAll = size(Phi)
	Ifi = nG+1; Ifip = Ifi+1; Ifim = Ifi-1
	Ila = ImaAll - nG; Ilap = Ila+1; Ilam = Ila-1
	Jfi = nG+1; Jfip = Jfi+1; Jfim = Jfi-1
	Jla = JmaAll - nG; Jlap = Jla+1; Jlam = Jla-1
    for k in 1:nPhi
        for j in Jfi:Jla
            @inbounds Phi[Ifim, j, k] = Phi[Ila, j, k]
            @inbounds Phi[Ilap, j, k] = Phi[Ifi, j, k]
        end
        for i in Ifi:Ila
            @inbounds Phi[i, Jfim, k] = Phi[i, Jla, k]
            @inbounds Phi[i, Jlap, k] = Phi[i, Jfi, k]
        end
    end
end