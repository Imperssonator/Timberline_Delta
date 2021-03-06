function out = Stackemup(DPdist)

%% Stackemup
% This function takes a degree of polymerization distribution as an input, in the form of an Nx2
% matrix:
% [DoP   # of chains]
% The output will be a histogram of the tie chain lengths


% as an initial example, we will use a distribution with average DP of 400
% (Mn = 67200) and a standard deviation of 50 (this is not the same as PDI)
% and assumes the chains follow a normal distribution, which is usually not
% the case

pi_length = 500;           % how many chains will we stack
DPindices = discretesample(DPdist(:,3),pi_length)';  % samples the DP distribution and returns a vector of the indices of which DPs are each new chain

Stack = initiate_stack(DPindices,DPdist);     % stack up the first two polymers

for i = 3:pi_length
    Stack = add_chain(DPdist(DPindices(i),1),Stack);
end

out = Stack;
stackplot(Stack)

end

function out = initiate_stack(DPindices,DPdist)

%% initiate stack
% takes the indices matrix, which specifies the locations of the DPs in the
% order in which they will stack, and takes the DPs from DPdist of those
% first two indices

DP1 = DPdist(DPindices(1),1);    % how long is each chain
DP2 = DPdist(DPindices(2),1);

hit1 = randi(DP1);          % which monomer of each chain collides with which
hit2 = randi(DP2);

out = [-hit1+1 DP1-hit1;...
       -hit2+1 DP2-hit2];
end

function out = add_chain(DP,Stack)

%% Add chain
% This function takes a given DP of the new chain and the existing stack in
% the form of
% [chain 1 start, chain 1 end
%  chain 2 start, chain 2 end...]
% and adds the new chain to either the beginning or end of the stack, with
% the collision happening on any monomer of the new chain and at any point
% on the overlapping region of the previous two chains in the existing
% stack

% the function returns the new stack with the added chain

side = randi(2);

if side == 2
    new_chain = collide(DP,Stack(end,:),Stack(end-1,:));
    out = [Stack; new_chain];
else
    new_chain = collide(DP,Stack(1,:),Stack(2,:));
    out = [new_chain; Stack];
end

end

function out = collide(DP,front,support)

%% collide
% takes a chain of length DP (integer scalar), and collides a random
% monomer with the frontier chain (front, 1x2 vector) and any location
% where front overlaps with support (1x2 vector)

% the output is a 1x2 vector

overlap = [max(front(1),support(1)), min(front(2),support(2))]; % find overlap of previous two chains

hitstack = randi(overlap(2)-overlap(1)+1)+overlap(1)-1; % where does new chain hit stack
hitchain = randi(DP);                                      % which monomer of new chain hits stack

out = [hitstack-hitchain+1, hitstack+DP-hitchain];

end