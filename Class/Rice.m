classdef Rice
    %RICE This class holds all the parameters for the Rician fading model.
    % It calculates theoretical envelope and phase PDFs
    % It does a Monte Carlo simulation using the parameters
    
    properties(Constant, Hidden = true)
        NumSamples = 2E6; % number of samples
        r = 0:0.001:6 % envelope range for PDF ploteter
        theta = -pi:0.001:pi; % phase range for PDF plot
    end 
    
    properties(Access = public)
        K; % Rician K factor
        r_hat_2; % root mean square of the signal
        phi; % phase parameter
    end
    
    properties(Hidden = true) 
        multipathFading; % Found based on the inputs
        envelopeProbability; % Calculated theoretical envelope probability
        phaseProbability; % Calculated theoretical phase probability
        xdataEnv; % Simulated envelope density plot x values 
        ydataEnv; % Simlated envelope density plot y valyes
        xdataPh; % Simulated phase density plot x values 
        ydataPh; % Simlated phase density plot y valyes
    end
    
    methods(Access = public)
        function obj = Rice(K,r_hat_2,phi)
            %ADDITIVESHADOWRICE Construct an instance of this class
            
            %   Assigning input values
            obj.K = input_Check(obj,K,'K',0,50);
            obj.r_hat_2 = input_Check(obj,r_hat_2,'\hat{r}^2',0.5,2.5);
            obj.phi = input_Check(obj,phi,'\phi',-pi,pi);
            
            % other calculated properties
            obj.multipathFading = complex_Multipath_Fading(obj);
            obj.envelopeProbability = envelope_PDF(obj);
            obj.phaseProbability = phase_PDF(obj);
            [obj.xdataEnv, obj.ydataEnv] = envelope_Density(obj);
            [obj.xdataPh, obj.ydataPh] = phase_Density(obj);
            
        end
    end
    
    methods(Access = private)
        
        function data = input_Check(obj, data, name, lower, upper) 
            % intput_Check checks the user inputs and throws errors
            
            % checks if input is empty
            if isempty(data)
                error(strcat(name,' must be a numeric input'));
            end
            
            % inputs must be a number
            if ~isnumeric(data)
               error(strcat(name,' must be a number, not a %s.', class(data)));
            end
            
            % input must be within the range
            if data < lower || data > upper
               error(strcat(name,' must be in the range [',num2str(lower),', ',num2str(upper),'].'));
            end
        end
        
        function [p, q] = means(obj)
            %means Calculates the means of the complex Gaussians 
            %representing the in-phase and quadrature components.

            p = sqrt(obj.K/(1+obj.K)) .* sqrt(obj.r_hat_2) .* cos(obj.phi);
            q = sqrt(obj.K/(1+obj.K)) .* sqrt(obj.r_hat_2) .* sin(obj.phi);

        end
        
        function [sigma] = scattered_Component(obj)
            %scattered_Component Calculates the power of the scattered 
            %signal component.    
            
            sigma = sqrt(obj.r_hat_2 ./(2 * (1 + obj.K)));
        
        end
        
        function [gaussians] = generate_Gaussians(obj, mean, sigma) 
            %generate_Gaussians Generates the Gaussian random variables 
            
            gaussians = normrnd(mean,sigma,[1,obj.NumSamples]);
        end
        
        function [multipathFading] = complex_Multipath_Fading(obj) 
            %complex_MultipathFading Generates the Rician fading model 
            
            [p, q] = means(obj);
            [sigma] = scattered_Component(obj);
            
            multipathFading = generate_Gaussians(obj, p, sigma) + 1i.* generate_Gaussians(obj, q, sigma);
        end    
        
        function [eProbTheor] = envelope_PDF(obj)
            %envelope_PDF Calculates the theoretical envelope PDF
            
            eProbTheor = 2 .* (1+obj.K) .* obj.r ./ obj.r_hat_2 ...
                .* exp(-obj.K - ((1+obj.K).*obj.r.^2 ./ obj.r_hat_2))...
                .* besseli(0, 2.*obj.r.*sqrt(obj.K.*(obj.K+1)./obj.r_hat_2));
        end
        
        function [pProbTheor] = phase_PDF(obj)
            %envelope_PDF Calculates the theoretical phase PDF
            
            pProbTheor = 1./(2*pi).*exp(-obj.K) .*(1 + (sqrt(4.*pi.*obj.K)...
                .* exp(obj.K.*(cos(obj.theta-obj.phi)).^2) .* cos(obj.theta-obj.phi))...
                .* (1 - qfunc(sqrt(2.*obj.K).*cos(obj.theta-obj.phi))));
        end
        
        function [xdataEnv, ydataEnv] = envelope_Density(obj)
            %envelope_Density Evaluates the envelope PDF
            R = sqrt((real(obj.multipathFading)).^(2) + (imag(obj.multipathFading)).^(2));

            [f,x] = ecdf(R);
            [ydataEnv, xdataEnv] = ecdfhist(f,x, 0:0.05:max(obj.r));
        end
        
        function [xdataPh, ydataPh] = phase_Density(obj)
            %envelope_Density Evaluates the envelope PDF
            T = angle((real(obj.multipathFading)) + (sqrt(-1)).*(imag(obj.multipathFading)));

            [f,x] = ecdf(T);
            [ydataPh, xdataPh] = ecdfhist(f,x, min(obj.theta)+0.03:0.06:max(obj.theta));
        end
            
    end
end

