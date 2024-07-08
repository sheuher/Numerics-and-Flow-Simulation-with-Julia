function log(Ima,Jma,nG,CFL,delta,DeltaT,tMax,nTMax,D,saveNT)
    println("""
    __________________________
    < PassiveScalarTransport2D >
    --------------------------

                ^__^
                (oo)_______
                (__)_       ) =3
                    ||----w |
                    ||     ||

    ---------------------------
    Ima$(lpad(Ima, 24, ' '))
    Jma$(lpad(Jma, 24, ' '))
    nG$(lpad(nG, 25, ' '))

    CFL$(lpad(CFL, 24, ' '))
    delta$(lpad(delta, 22, ' '))
    DeltaT$(lpad(DeltaT, 21, ' '))
    tMax$(lpad(tMax, 23, ' '))
    nTMax$(lpad(nTMax, 22, ' '))
    D$(lpad(D, 26, ' '))
    saveNT$(lpad(saveNT, 21, ' '))
    ---------------------------

    """)
end

function log(I,DeltaT,t)
    println("|$(rpad(I, 6, ' '))|$(rpad(DeltaT, 6, ' '))$(lpad("|$t", 10, ' '))|")
end

function log()
    println("---------------------------")
end