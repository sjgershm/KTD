function results = KTD_sim(sim,param)
    
    results.sim = sim;
    if nargin < 2; param = KTD_defparam; end
    
    switch sim
        
        case 1
            nTrials = 10;
            trial_length = 10;
            dur = 4;
            ons = 3;
            basis = 'CSC';
            
            % x = [X,Z,A,B,Y]
            
            % AX+
            x1 = zeros(trial_length,5);
            x1(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x1(:,3) = KTD_make_stimulus(ons,dur,trial_length);
            f1 = construct_basis(basis,x1);
            r1 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % BY+
            x2 = zeros(trial_length,5);
            x2(:,4) = KTD_make_stimulus(ons,dur,trial_length);
            x2(:,5) = KTD_make_stimulus(ons,dur,trial_length);
            f2 = construct_basis(basis,x2);
            r2 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % A-
            x3 = zeros(trial_length,5);
            x3(:,3) = KTD_make_stimulus(ons,dur,trial_length);
            f3 = construct_basis(basis,x3);
            r3 = zeros(trial_length,1);
            
            % Z->X
            x4 = zeros(trial_length,5);
            x4(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x4(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f4 = construct_basis(basis,x4);
            r4 = zeros(trial_length,1);
            
            F = [repmat(f1,nTrials,1); repmat(f2,nTrials,1); repmat(f3,nTrials,1);  repmat(f4,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r2,nTrials,1); repmat(r3,nTrials,1);  repmat(r4,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model)
                w(n,:)=model(n).w;
                c(n,:)=model(n).C(ons,:);
                k(n,:)=model(n).K;
            end
            results.model = model;
            results.W = w(1:trial_length:end,trial_length+1);
            results.C = c(3:trial_length:end,:);
            results.K = k(3:trial_length:end,:);
            
        case 2
            nTrials = 10;
            trial_length = 10;
            dur = 4;
            ons = 3;
            basis = 'CSC';
            
            % x = [X,Z,A,B,Y]
            
            % AX+
            x1 = zeros(trial_length,5);
            x1(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x1(:,3) = KTD_make_stimulus(ons,dur,trial_length);
            f1 = construct_basis(basis,x1);
            r1 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % BY+
            x2 = zeros(trial_length,5);
            x2(:,4) = KTD_make_stimulus(ons,dur,trial_length);
            x2(:,5) = KTD_make_stimulus(ons,dur,trial_length);
            f2 = construct_basis(basis,x2);
            r2 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % B-
            x3 = zeros(trial_length,5);
            x3(:,4) = KTD_make_stimulus(ons,dur,trial_length);
            f3 = construct_basis(basis,x3);
            r3 = zeros(trial_length,1);
            
            % Z->X
            x4 = zeros(trial_length,5);
            x4(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x4(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f4 = construct_basis(basis,x4);
            r4 = zeros(trial_length,1);
            
            F = [repmat(f1,nTrials,1); repmat(f2,nTrials,1); repmat(f3,nTrials,1);  repmat(f4,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r2,nTrials,1); repmat(r3,nTrials,1);  repmat(r4,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model)
                w(n,:)=model(n).w;
                c(n,:)=model(n).C(ons,:);
                k(n,:)=model(n).K;
            end
            results.model = model;
            results.W = w(1:trial_length:end,trial_length+1);
            results.C = c(3:trial_length:end,:);
            results.K = k(3:trial_length:end,:);
            
        case 3
            nTrials = 10;
            trial_length = 10;
            dur = 4;
            ons = 3;
            basis = 'CSC';
            
            % x = [A,B]
            
            % BA+
            x1 = zeros(trial_length,2);
            x1(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x1(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f1 = construct_basis(basis,x1);
            r1 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % A-
            x2 = zeros(trial_length,2);
            x2(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            f2 = construct_basis(basis,x2);
            r2 = zeros(trial_length,1);
            
            % null
            x3 = zeros(trial_length,2);
            f3 = construct_basis(basis,x3);
            r3 = zeros(trial_length,1);
            
            F = [repmat(f1,nTrials,1); repmat(f2,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r2,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(1:trial_length:end,trial_length+1);
            results(1).model = model;
            results(1).W = W;
            
            F = [repmat(f1,nTrials,1); repmat(f3,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r3,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(1:trial_length:end,trial_length+1);
            results(2).model = model;
            results(2).W = W;
            
        case 4
            
            nTrials = 10;
            trial_length = 10;
            dur = 4;
            ons = 3;
            basis = 'CSC';
            
            % x = [A,B]
            
            % BA+
            x1 = zeros(trial_length,2);
            x1(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x1(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f1 = construct_basis(basis,x1);
            r1 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % A-
            x2 = zeros(trial_length,2);
            x2(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            f2 = construct_basis(basis,x2);
            r2 = zeros(trial_length,1);
            
            % null
            x3 = zeros(trial_length,2);
            f3 = construct_basis(basis,x3);
            r3 = zeros(trial_length,1);
            
            F = [repmat(f2,nTrials,1); repmat(f1,nTrials,1)];
            r = [repmat(r2,nTrials,1); repmat(r1,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(3:trial_length:end,trial_length+1);
            results(1).model = model;
            results(1).W = W;
            
            F = [repmat(f3,nTrials,1); repmat(f1,nTrials,1)];
            r = [repmat(r3,nTrials,1); repmat(r1,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(3:trial_length:end,trial_length+1);
            results(2).model = model;
            results(2).W = W;
            
        case 5
            % serial compound latent inhibition
            
            s0 = param.s;
            s = [s0 s0 s0*2 s0*2];
            lr0 = param.lr;
            lr = [lr0 lr0 lr0/2 lr0/2];
            n2 = 10;
            n1 = 10;
            buffer = 5;
            X{1} = [zeros(n1,1) ones(n1,1); zeros(buffer,2); ones(n2,1) zeros(n2,1)];
            X{2} = [ones(n1,1) zeros(n1,1); zeros(buffer,2); ones(n2,1) zeros(n2,1)];
            X{3} = X{2};
            X{4} = X{2};
            r = [zeros(n1+buffer,1); ones(n2,1)];
            for i = 1:4
                if i < 4
                    param.lr = [zeros(n1+buffer,1)+lr(i); zeros(n2,1)+lr0];
                    param.s = [zeros(n1+buffer,1)+s(i); zeros(n2,1)+s0];
                else
                    param.lr = [zeros(n1+buffer,1)+lr0; zeros(n2,1)+lr(i)];
                    param.s = [zeros(n1+buffer,1)+s0; zeros(n2,1)+s(i)];
                end
                model = kalmanTD(X{i},r,param);
                for n=1:length(model); w(n,1)=model(n).w(1); end
                results.W(:,i) = w(end-9:end);
            end
            
        case 6
            % latent inhibition
            
            n = 10;
            X{1} = [ones(n,1); ones(n,1)];
            X{2} = [zeros(n,1); ones(n,1)];
            r = [zeros(n,1); ones(n,1)];
            
            for i = 1:length(X)
                model = kalmanRW(X{i},r);
                results.W(:,i) = [model.w0];
                results.k(:,i) = [model.K];
            end
            
        case 7
            % overshadowing
            n = 10;
            X{1} = [ones(n,2); zeros(n,2)];
            X{2} = [ones(n,2); ones(n,1) zeros(n,1)];
            r = [ones(n,1); zeros(n,1)];
            
            for i = 1:length(X)
                model = kalmanRW(X{i},r);
                results.W(i,:) = model(end).w0;
                results.k(i,:) = model(end).K;
            end
            
        case 8
            % forward blocking
            n = 10;
            X{1} = [ones(n,1) zeros(n,1); ones(n,2); zeros(n,2)];
            X{2} = [ones(n,1) zeros(n,1); ones(n,2); ones(n,1) zeros(n,1)];
            r = [ones(n*2,1); zeros(n,1)];
            
            for i = 1:length(X)
                model = kalmanRW(X{i},r);
                results.W(i,:) = model(end).w0;
                results.k(i,:) = model(end).K;
            end
            
        case 9
            % overexpectation
            n = 10;
            X{1} = [ones(n,1) zeros(n,1); zeros(n,1) ones(n,1); ones(n,2); zeros(n,2)];
            X{2} = [ones(n,1) zeros(n,1); zeros(n,1) ones(n,1); ones(n,2); ones(n,1) zeros(n,1)];
            r = [ones(n*3,1); zeros(n,1)];
            
            for i = 1:length(X)
                model = kalmanRW(X{i},r);
                results.W(i,:) = model(end).w0;
                results.k(i,:) = model(end).K;
            end
            
        case 10
            % conditioned inhibition
            n = 12;
            X{1} = [repmat([1 0; 1 1],n,1); zeros(n,2)];
            X{2} = [repmat([1 0; 1 1],n,1); ones(n,1) zeros(n,1)];
            r = [repmat([1; 0; 1; 0; 1; 0; 1; 1],n/4,1); zeros(n,1)];
            
            for i = 1:length(X)
                model = kalmanRW(X{i},r);
                results.W(i,:) = model(end).w0;
                results.k(i,:) = model(end).K;
            end
            
        case 11
            % serial compound backward blocking
            nTrials = 10;
            trial_length = 10;
            dur = 4;
            ons = 3;
            basis = 'CSC';
            
            % x = [A,B]
            
            % BA+
            x1 = zeros(trial_length,2);
            x1(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            x1(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f1 = construct_basis(basis,x1);
            r1 = KTD_make_stimulus(ons+dur-1,1,trial_length);
            
            % A+
            x2 = zeros(trial_length,2);
            x2(:,1) = KTD_make_stimulus(ons,dur,trial_length);
            f2 = construct_basis(basis,x2);
            r2 = r1;
            
            % B+
            x3 = zeros(trial_length,2);
            x3(:,2) = KTD_make_stimulus(1,dur,trial_length);
            f3 = construct_basis(basis,x3);
            r3 = KTD_make_stimulus(dur,1,trial_length);
            
            F = [repmat(f1,nTrials,1); repmat(f2,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r2,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(1:trial_length:end,trial_length+1); % wB
            results(1).model = model;
            results(1).W = W;
            
            F = [repmat(f1,nTrials,1); repmat(f3,nTrials,1)];
            r = [repmat(r1,nTrials,1); repmat(r3,nTrials,1)];
            
            model = kalmanTD(F,r,param);
            for n=1:length(model); w(n,:)=model(n).w; end
            W = w(ons:trial_length:end,3);    % wA
            results(2).model = model;
            results(2).W = W;
            
        case 12
            % CS-US interval, fixed ITI
            
            nTrials = 20;
            trial_length = 25;
            dur = 1:4;
            ons = 1;
            basis = 'CSC';
            
            for i = 1:length(dur)
                x1 = KTD_make_stimulus(ons,dur(i),trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(ons+dur(i)-1,1,trial_length);
                
                F = repmat(f1,nTrials,1);
                r = repmat(r1,nTrials,1);
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length:end,1); % wB
                clear w
                results(i).model = model;
                results(i).W = W;
            end
            
        case 13
            % CS-US interval, fixed ISI/ITI ratio
            
            nTrials = 20;
            %dur = 1:4;
            dur = [2 2 2 2];
            %trial_length = 10*dur;
            trial_length = [10 20 30 40];
            ons = 1;
            basis = 'CSC';
            
            for i = 1:length(dur)
                x1 = KTD_make_stimulus(ons,dur(i),trial_length(i));
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(ons+dur(i)-1,1,trial_length(i));
                
                F = repmat(f1,nTrials,1);
                r = repmat(r1,nTrials,1);
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length(i):end,1); % wB
                clear w
                results(i).model = model;
                results(i).W = W;
            end
            
        case 14
            % CS duration effect
            
            nTrials = 100;
            trial_length = 10;
            dur = [2 4 6];
            ons = 1;
            basis = 'CSC';
            
            % x = [A,B]
            
            for i = 1:length(dur)
                
                % AB+
                x1 = zeros(trial_length,2);
                x1(:,1) = KTD_make_stimulus(ons,dur(i),trial_length);
                x1(:,2) = KTD_make_stimulus(ons,dur(i),trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(ons+dur(i)-1,1,trial_length);
                
                % B-
                x2 = zeros(trial_length,2);
                x2(:,2) = KTD_make_stimulus(1,dur(i),trial_length);
                f2 = construct_basis(basis,x2);
                r2 = r1;
                
                F = [repmat(f1,nTrials,1); repmat(f2,nTrials,1)];
                r = [repmat(r1,nTrials,1); repmat(r2,nTrials,1)];
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length:end,1);
                clear w
                results(i).model = model;
                results(i).W = W;
            end
            
        case 15
            % second-order simultaneous vs. serial conditioning
            
            nTrials = 3;
            rep = [1 2 3];
            trial_length = 10;
            dur = 4;
            onsA = 3;
            onsB = [onsA-1 onsA];
            basis = 'CSC';
            param.TD = 0;
            
            for i = 1:length(rep)
                for j = 1:length(onsB)
                    
                    % x = [A,B]
                    
                    % A+
                    x1 = zeros(trial_length,2);
                    x1(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                    f1 = construct_basis(basis,x1);
                    r1 = KTD_make_stimulus(onsA+dur-1,1,trial_length);
                    
                    % B->A-
                    x2 = zeros(trial_length,2);
                    x2(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                    x2(:,2) = KTD_make_stimulus(onsB(j),dur,trial_length);
                    f2 = construct_basis(basis,x2);
                    r2 = zeros(trial_length,1);
                    
                    f2 = repmat(f2,rep(i),1);
                    r2 = repmat(r2,rep(i),1);
                    
                    F = repmat([f1; f2],nTrials,1);
                    r = repmat([r1; r2],nTrials,1);
                    
                    model = kalmanTD(F,r,param);
                    for n=1:length(model); w(n,:)=model(n).w; end
                    W = w(1:trial_length:end,trial_length+onsB(j));
                    results(i,j).model = model;
                    results(i,j).W = W;
                    results(i,j).rep = rep;
                end
            end
            
        case 16
            
            % second-order latent inhibition
            
            nTrials = 4;
            trial_length = 10;
            dur = 4;
            onsA = 2;
            onsB = 2;
            basis = 'CSC';
            pre = 8;
            
            for j = 1:2
                
                % x = [A,B]
                
                % B-
                x0 = zeros(trial_length,2);
                if j == 1
                    x0(:,2) = KTD_make_stimulus(onsB,dur,trial_length);
                end
                f0 = construct_basis(basis,x0);
                r0 = zeros(trial_length,1);
                
                % A+
                x1 = zeros(trial_length,2);
                x1(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(onsA+dur-1,1,trial_length);
                f1 = repmat(f1,2,1);
                r1 = repmat(r1,2,1);
                
                % BA-
                x2 = zeros(trial_length,2);
                x2(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                x2(:,2) = KTD_make_stimulus(onsB,dur,trial_length);
                f2 = construct_basis(basis,x2);
                r2 = zeros(trial_length,1);
                
                F = [repmat(f0,pre,1); repmat([f1; f2],nTrials,1)];
                r = [repmat(r0,pre,1); repmat([r1; r2],nTrials,1)];
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length:end,trial_length+onsB);
                results(j).model = model;
                results(j).W = W;
            end
            
        case 17
            
            % second-order extinction
            
            nTrials = 2;
            trial_length = 10;
            dur = 2;
            onsA = 3;
            onsB = 1;
            basis = 'CSC';
            ext = 8;
            
            for j = 1:2
                
                % x = [A,B]
                
                % A+
                x1 = zeros(trial_length,2);
                x1(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(onsA+dur-1,1,trial_length);
                f1 = repmat(f1,2,1);
                r1 = repmat(r1,2,1);
                
                % B->A-
                x2 = zeros(trial_length,2);
                x2(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                x2(:,2) = KTD_make_stimulus(onsB,dur,trial_length);
                f2 = construct_basis(basis,x2);
                r2 = zeros(trial_length,1);
                
                % A-
                x3 = zeros(trial_length,2);
                if j == 1
                    x3(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                end
                f3 = construct_basis(basis,x3);
                r3 = zeros(trial_length,1);
                
                F = [repmat([f1; f2],nTrials,1); repmat(f3,ext,1)];
                r = [repmat([r1; r2],nTrials,1); repmat(r3,ext,1)];
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length:end,trial_length+onsB);
                results(j).model = model;
                results(j).W = W;
            end
            
        case 18
            
            % serial compound overshadowing
            
            nTrials = 2;
            trial_length = 10;
            dur = 2;
            onsA = 3;
            onsB = 1;
            basis = 'CSC';
            ext = 8;
            
            for j = 1:2
                
                % x = [A,B]
                
                % A+
                x1 = zeros(trial_length,2);
                x1(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(onsA+dur-1,1,trial_length);
                f1 = repmat(f1,2,1);
                r1 = repmat(r1,2,1);
                
                % B->A-
                x2 = zeros(trial_length,2);
                x2(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                x2(:,2) = KTD_make_stimulus(onsB,dur,trial_length);
                f2 = construct_basis(basis,x2);
                r2 = zeros(trial_length,1);
                
                % A-
                x3 = zeros(trial_length,2);
                if j == 1
                    x3(:,1) = KTD_make_stimulus(onsA,dur,trial_length);
                end
                f3 = construct_basis(basis,x3);
                r3 = zeros(trial_length,1);
                
                F = [repmat([f1; f2],nTrials,1); repmat(f3,ext,1)];
                r = [repmat([r1; r2],nTrials,1); repmat(r3,ext,1)];
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(1:trial_length:end,trial_length+onsB);
                results(j).model = model;
                results(j).W = W;
            end
            
        case 19
            % mediated extinction and recovery from overshadowing
            
            nTrials = 10;
            trial_length = 10;
            dur = 3;
            ons2 = [1 3];
            ons1 = 1;
            basis = 'CSC';
            
            for j = 1:length(ons2)
                
                % x = [A,X,B,Y]
                
                % AX+
                x1 = zeros(trial_length,4);
                x1(:,1) = KTD_make_stimulus(ons1,dur,trial_length);
                x1(:,2) = KTD_make_stimulus(ons2(j),dur,trial_length);
                f1 = construct_basis(basis,x1);
                r1 = KTD_make_stimulus(ons2(j)+dur-1,1,trial_length);
                
                % BY+
                x2 = zeros(trial_length,4);
                x2(:,3) = KTD_make_stimulus(ons1,dur,trial_length);
                x2(:,4) = KTD_make_stimulus(ons2(j),dur,trial_length);
                f2 = construct_basis(basis,x2);
                r2 = KTD_make_stimulus(ons2(j)+dur-1,1,trial_length);
                
                % A-
                x3 = zeros(trial_length,4);
                x3(:,1) = KTD_make_stimulus(ons1,dur,trial_length);
                f3 = construct_basis(basis,x3);
                r3 = zeros(trial_length,1);
                
                % X-
                x4 = zeros(trial_length,4);
                x4(:,2) = KTD_make_stimulus(ons2(j),dur,trial_length);
                f4 = construct_basis(basis,x4);
                r4 = -KTD_make_stimulus(ons2(j)+dur-1,1,trial_length);
                
                % Y-
                x5 = zeros(trial_length,4);
                x5(:,4) = KTD_make_stimulus(ons2(j),dur,trial_length);
                f5 = construct_basis(basis,x5);
                r5 = -KTD_make_stimulus(ons2(j)+dur-1,1,trial_length);
                
                F = [repmat([f1; f2],nTrials,1); repmat(f3,nTrials,1)];
                r = [repmat([r1; r2],nTrials,1); repmat(r3,nTrials,1)];
                
                model = kalmanTD(F,r,param);
                for n=1:length(model); w(n,:)=model(n).w; end
                W = w(ons2(j):trial_length:end,[trial_length+ons2(j) 3*trial_length+ons2(j)]);
                clear w
                results(j).model = model;
                results(j).W = W;
            end
            
    end
    