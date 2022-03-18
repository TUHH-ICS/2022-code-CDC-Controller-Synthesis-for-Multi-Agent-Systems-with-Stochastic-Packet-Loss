classdef DisturbedAgent < LinearFrictionMass
    %DISTURBEDAGENT Examplary agent implementation with simple LTI
    %dynamics that performs a formation control manoeuvre.
    
    properties(GetAccess = public, SetAccess = private)
        disturbance % Disturbance signal of the agents
    end
    
    properties(GetAccess = public, SetAccess = immutable)
        ref
        neighbours
    end
    
    properties(GetAccess = private, SetAccess = immutable)
        controller
    end
    
    methods
        function obj = DisturbedAgent(id, initialPos, ref, neighbours, disturbed)
            %DISTURBEDAGENT Construct an instance of this class
            %   Sets up the correct agent dynamics and initializes the
            %   agent and consensus protocol to the given initial position.
            
            data = load('controller');
            
            initialVel = zeros(size(initialPos));
            obj@LinearFrictionMass(id, data.dT, data.m, data.b, initialPos, initialVel);
            obj.ref = ref;
            obj.neighbours = neighbours;
            
            % Load controller matrices. We cannot have Bc_K & Dc_K, since
            % the plant is distributed at the output
            [Ad_K, Bd_K, Cd_K, Dd_K] = ssdata(data.Kd);
            [Ac_K, ~,    Cc_K, ~]    = ssdata(data.Kc);
            
            % Initialize controller
            obj.controller = DiscreteLtiDynamics(Ad_K, [Bd_K, Ac_K],...
                                                 Cd_K, [Dd_K, Cc_K]);
            
            % Set disturbance status
            obj.disturbance = zeros(obj.dim, 1);
            if disturbed
                obj.disturbance(disturbed) = 1;
            end
        end
        
        function step(obj)
            % Receive messages from the network
            messages = obj.receive();
            
            % Filter to only neighbours
            if ~isempty(messages)
                mask = any(obj.neighbours == [messages.sender]', 2);
                messages = messages(mask);
            end
            
            if ~isempty(messages)               
                % Calculate new formation reference. We use the standard
                % Laplacian, therefore we sum up the data
                data = [messages.data];
                err  = sum(obj.position - [data.position] - obj.ref, 2);
                ctr  = sum(obj.controller.x - [data.ctr_state], 2);
            else
                err = zeros(obj.dim,1);
                ctr = zeros(size(obj.controller.A,1),1);
            end
            
            % Evaluate controller & agent dynamics
            u = obj.controller.step([-err; ctr]);
            obj.move(u + obj.disturbance);
            
            % Apply disturbance only once
            obj.disturbance(:) = 0;
            
            % Send message to network, include only the position 
            data = struct;
            data.position  = obj.position - obj.ref;
            data.ctr_state = obj.controller.x;
            obj.send(data)
        end
    end
end