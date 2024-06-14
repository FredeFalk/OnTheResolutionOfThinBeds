function [mods] = write_modfile(N,RES1,RES2,THICK1,THICK2,Fmode,wedge_t,wedge_r)
mods = zeros(N,4);
fileID = fopen('tempo.mod','w');

%Setting up file head
fprintf(fileID,'%6s\n','Model filer, Forward SkyTEM (free text line)');
fprintf(fileID,'%3.0f %3.0f\n',N,0);


for i = 1:N
    fprintf(fileID,'%3.0f %3.0f %18s\n',i,1,'SkyTEM304_combined.tem');
end

fprintf(fileID,'%0.0f\n',-1);

%add N random models
if Fmode == 0
    for i = 1:N

        R1 = power(10,drawrandom(RES1));
        R2 = power(10,drawrandom(RES2));
        R3 = R1;
        TK1 = drawrandom(THICK1);
        TK2 = drawrandom(THICK2);



        model = [R1,R2,R3,TK1,TK2];
        mods(i,:) = [R1 R2 TK1 TK2];

        fprintf(fileID,'%0.0f\n',3);
        fprintf(fileID,'%9.2f %4.0f\n',model(1),-1);  %RES1
        fprintf(fileID,'%9.2f %4.0f\n',model(2),-1);   %RES2
        fprintf(fileID,'%9.2f %4.0f\n',model(3),-1);  %RES3
        fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %THK1
        fprintf(fileID,'%9.2f %4.0f\n',model(5),-1);   %THK2
        fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %DEP1
        fprintf(fileID,'%9.2f %4.0f\n',model(5)+model(4),-1);  %DEP2

    end
else
    if Fmode == 1
        display('Gridding On')
        for i = 1:sqrt(N)
            for j = 1:sqrt(N)

                %Define Grid
                R1 = power(10,RES1(1)+(RES1(2)-RES1(1))*((i-1)/(sqrt(N)-1)));
                R2 = power(10,RES2(1)+(RES2(2)-RES2(1))*((i-1)/(sqrt(N)-1)));
                R3 = R1;

                TK1 = THICK1(1)+(THICK1(2)-THICK1(1))*((i-1)/(sqrt(N)-1));
                TK2 = THICK2(1)+(THICK2(2)-THICK2(1))*((j-1)/(sqrt(N)-1));


                model = [R1,R2,R3,TK1,TK2];

                mods(j+(i-1)*sqrt(N),:) = [R1 R2 TK1 TK2];

                fprintf(fileID,'%0.0f\n',3);
                fprintf(fileID,'%9.2f %4.0f\n',model(1),-1);  %RES1
                fprintf(fileID,'%9.2f %4.0f\n',model(2),-1);   %RES2
                fprintf(fileID,'%9.2f %4.0f\n',model(3),-1);  %RES3
                fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %THK1
                fprintf(fileID,'%9.2f %4.0f\n',model(5),-1);   %THK2
                fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %DEP1
                fprintf(fileID,'%9.2f %4.0f\n',model(5)+model(4),-1);  %DEP2
            end
        end
    else
        if Fmode == 2
            for i = 1:sqrt(N)
                for j = 1:sqrt(N)

                    R1 = power(10,RES1(i));
                    R2 = power(10,RES2(i));
                    R3 = R1;
                    TK1 = THICK1(i);
                    TK2 = THICK2(j);

                    model = [R1,R2,R3,TK1,TK2];

                    mods((i-1)*sqrt(N)+j,:) = [R1 R2 TK1 TK2];

                    fprintf(fileID,'%0.0f\n',3);
                    fprintf(fileID,'%9.2f %4.0f\n',model(1),-1);  %RES1
                    fprintf(fileID,'%9.2f %4.0f\n',model(2),-1);   %RES2
                    fprintf(fileID,'%9.2f %4.0f\n',model(3),-1);  %RES3
                    fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %THK1
                    fprintf(fileID,'%9.2f %4.0f\n',model(5),-1);   %THK2
                    fprintf(fileID,'%9.2f %4.0f\n',model(4),-1);   %DEP1
                    fprintf(fileID,'%9.2f %4.0f\n',model(5)+model(4),-1);  %DEP2
                end
            end
        end
    end

    if Fmode == 3
        for i = 1:N/size(wedge_t,2)

            Resi = wedge_r;

            for j = 1:size(wedge_t,2)

                Thicki = wedge_t(2*i-1:2*i,j);

                if Thicki(2) < 0.02
                    fprintf(fileID,'%0.0f\n',2);
                    fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                    fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);   %RES3

                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
                else

                    fprintf(fileID,'%0.0f\n',3);
                    fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                    fprintf(fileID,'%9.2f %4.0f\n',Resi(2),-1);   %RES2
                    fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);  %RES3

                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(2),-1);   %THK2

                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
                    fprintf(fileID,'%9.2f %4.0f\n',Thicki(2)+Thicki(1),-1);  %DEP2
                end
            end
        end
    end

    if Fmode == 4
        for i = 1:N

            Resi = wedge_r(i,:);
            Thicki = wedge_t(i,:);
            mods(i,:) = [Resi(1) Resi(2) Thicki(1) Thicki(2)];

            if Thicki(2) < 0.02
                fprintf(fileID,'%0.0f\n',2);
                fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);   %RES3

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
            else

                fprintf(fileID,'%0.0f\n',3);
                fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                fprintf(fileID,'%9.2f %4.0f\n',Resi(2),-1);   %RES2
                fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);  %RES3

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(2),-1);   %THK2

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(2)+Thicki(1),-1);  %DEP2
            end
        end
    end

    if Fmode == 100
        for i = 1:N

            Resi = wedge_r(i,:);
            Thicki = wedge_t(i,:);

            mods(i,:) = [Resi(1) Resi(2) Thicki(1) Thicki(2)];

            if Thicki(2) < 0.02
                fprintf(fileID,'%0.0f\n',2);
                fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);   %RES3

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
            else

                fprintf(fileID,'%0.0f\n',3);
                fprintf(fileID,'%9.2f %4.0f\n',Resi(1),-1);  %RES1
                fprintf(fileID,'%9.2f %4.0f\n',Resi(2),-1);   %RES2
                fprintf(fileID,'%9.2f %4.0f\n',Resi(3),-1);  %RES3

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %THK1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(2),-1);   %THK2

                fprintf(fileID,'%9.2f %4.0f\n',Thicki(1),-1);   %DEP1
                fprintf(fileID,'%9.2f %4.0f\n',Thicki(2)+Thicki(1),-1);  %DEP2
            end
        end
    end

    fclose(fileID);
end
end