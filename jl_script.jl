using Pkg
Pkg.add("QPDAS")
using QPDAS
using LinearAlgebra


# Funktion read_datei bekommt Weg zum Datei im System
# Datei in solchen Format angegeben
#    input: path to the file with rows {xi1}  {xi2} {yi}
# Funktion liest Datei, und gibt Antwort im Form
#    output:Array([xi1,xi2],yi)

function read_datei(filename)
    x = []
    y = []
    a = readlines(filename)

    for k=1:length(a)
        temp_row = split(a[k])
        x_vector = [parse(Float64, temp_row[1]),parse(Float64, temp_row[2])]
        push!(x,x_vector)
        push!(y,parse(Float64, temp_row[3]))
    end
    return (x,y)
end

# Berechnung vom Skalarprodukt von 2 Elementigen Vektoren
# Es wird bei der Matrix Berechnung benutzt
#   input: x1, 2-Vektor
#          x2, 2-Vektor
#    output: double, Ergebnis von skalarprodukt
function skalar_produkt(x1,x2)
    return x1[1]*x2[1]+x1[2]*x2[2]
end

# Berechnung von Matrix M = yi*yj*<xi*xj>, die im Qudratisches Programm benutzt wird
# input: Array([xi1,xi2],yi), aus read_datei
# output: k*k Matrix, wo k Anzahl von Biespielen
function matrix_berechnung(array)
    matrix = zeros(length(array[1]),length(array[1]))
    for i=1:length(array[1])
        for j = 1:length(array[1])
            matrix[i,j] = array[2][i]*array[2][j]*skalar_produkt(array[1][i],array[1][j])
        end
    end
    return matrix
end

### Creating a variables for Quadratic program
# P = M
# z = 1, vector
# b = 0 
# c = -I, Einheitsmatrix
# d = 0, vector
# A Matrix bei Diagonalelementen hat spalte y aus Beispiel

# calculate_write_w_b - Sammlung von obigen Funktionen, um Datei zu lesen und Werte zu speicher
#                       + bereitstellung von Quadratischen Programm und Loesung fuer w und b
#                       + schreiben von w und b als Ausgabe im Datei mit dem Name "Beispiel{i}_vektor.txt"
function calculate_write_w_b(filename)
    content = read_datei(filename)
    M = matrix_berechnung(content)
    P_my = M
    #print("Matrix M",M)
    z_my = zeros(Int64(sqrt(length(M))))
    for k=1:Int64(sqrt(length(M)))
        z_my[k] = Float64(-1)
    end

    b = zeros(1)#Int64(sqrt(length(M))))
    #for k=1:Int64(sqrt(length(M)))
    #    b[k] = Float64(0)
    #end
    #print("Vektor b",b)
    C = zeros(Int64(sqrt(length(M))),Int64(sqrt(length(M))))
    for i=1:Int64(sqrt(length(M)))
        C[i,i] = -1
    end
    #print("Matrix C",C)
    A = zeros(Int64(sqrt(length(M))))#,Int64(sqrt(length(M))))
    for k=1:Int64(sqrt(length(M)))
        A[k] = Float64(content[2][k])
    end
    A = transpose(A)

    A = Matrix{Float64}(A)
    #for i=1:Int64(sqrt(length(M)))
    #    for j = 1:Int64(sqrt(length(M)))
    #        if(i == j)
    #            A[i,j] = content[2][i]
    #        end
    #    end
    #end
    #print("Matrix A",A)
    d = zeros(Int64(sqrt(length(M))))
    

    P_my  
    
    P_my = P_my+1/1000*I 

    qp = QuadraticProgram(A, b, C, d, z_my, P_my)
    sol,val = solve!(qp)
    #print(sol)

    w = zeros(Int64(2))
    for i=1:(Int64(sqrt(length(M))))
        w += sol[i]*content[2][i]*content[1][i] 
    end


    
    I_len = Int64(sqrt(length(M)))
    b = 0
    for i=1:Int64(sqrt(length(M)))
        if(sol[i]<10^(-5))
            I_len -=1
            continue
        end
        temp = 0
        for j =1:Int64(sqrt(length(M)))
            if(sol[j]<10^(-5))
                continue
            end
            temp += sol[j]*content[2][j]*skalar_produkt(content[1][j],content[1][i]) 
        end
        b -= content[2][i] - temp
    end
    b = b/I_len

    datei = open(filename[1:length(filename)-4]*"_vektor.txt","w");
    str = string(w[1])*"\n"*string(w[2])*"\n"*string(b)
    write(datei,str);
    close(datei)
    return(w,b)
end