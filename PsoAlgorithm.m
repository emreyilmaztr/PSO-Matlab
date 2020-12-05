clear all;
tol=10.^(-6);
% self confidence. [1.5 - 2]
C1 = 1.7;
% swarm confidence. [2 - 2.5]
C2 = 2.2;
% intertia factor [0.4 - 1.4]
W = 0.7;
ParticleCount = 15;
MaxIteration = 500;
% J = particle number.
% i = iteration number.

f=@(x)(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2; %%Booth Function
% f=@(x)(1.5-x(1)*(1-x(2)))^2+(2.25-x(1)*(1-x(2)^2))^2+...
%     (2.625-x(1)*(1-x(2)^3))^2; %% Beale Function
n=2; %% dimension of the problem
% Initial V.
%  Vji = rand(1,ParticleCount);
  Vji = rand(ParticleCount,n);
% Initial X.
%  Xji = rand(1,ParticleCount);
  Xji = rand(ParticleCount,n);

% Personal Best.
for i=1:ParticleCount
PBest(i,:) = Xji(i,:);
result(i)=f(PBest(i,:));
end
[MinVal,MinIndex] = min(result);
% Group Best.

GBest = PBest(MinIndex,:);
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempGbest = GBest - 1;
GBestCount = 0;


iter = 1;
termState = 1;

while ( (GBestCount ~= 50) && (iter < MaxIteration) )
   
for i=1:ParticleCount

    Vji(i,:) = W * Vji(i,:) + C1 * rand(1,n).*(PBest(i,:)-Xji(i,:)) + C2 * rand(1,n) .*(GBest - Xji(i,:));
   
    Xji(i,:) = Xji(i,:) + Vji(i,:);
end
 
    for j = 1:ParticleCount
       
        oldValue = f(PBest(j,:));
        newValue = f(Xji(j,:));
        if (newValue < oldValue)
            PBest(j,:) = Xji(j,:);
        end    
    end
    for i=1:ParticleCount
         result(i)=f(PBest(i,:));
    end
    [MinVal,MinIndex] = min(result)
   GBest = PBest(MinIndex,:)
   
    if (TempGbest == GBest)
        GBestCount = GBestCount + 1;
    else
        GBestCount = 0;
    end
   
   TempGbest = GBest;
   G(iter,:)=GBest;
   S(iter)=f(GBest);
   if iter>=2
       h=abs(S(iter)-S(iter-1))
   else
       h=1;
   end
   if iter>40 && h<tol
       break
   end
%    plot(iter,S(iter))
%    hold on
    iter = iter + 1;
   
end
a=length(S)
Sonuc=S
%% Graphical View
% [y,t]=meshgrid(-4:0.1:4, -4:0.1:4);
% z=(y+2.*t-7).^2+(2.*y+t-5).^2; %% Booth Function
% % z=(1.5-y.*(1-t)).^2+(2.25-y.*(1-t.^2)).^2+(2.625-y.*(1-t.^3)).^2;
% meshc(y,t,z)
% hold on
% for i=1:a
% plot3(G(i,1),G(i,2),S(i),'.r','MarkerSize',15);
% hold on
% end
% title('Objective Function Values'); xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis')
%% Convergence Analysis
for i=1:length(S)
    x(i)=i;
end
plot(x,S,'LineWidth',2); grid on;
title('Objective Function Values'); xlabel('Iteration'); ylabel('f(x)');