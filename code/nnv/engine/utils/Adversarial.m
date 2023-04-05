classdef Adversarial
    %Adversarial class
    %   Contains static methods for adversarial robustness testing
    
    properties
    end
    
    methods(Static)
        
        %Adversarial testing for specific input vector on specific attacked
        %positions
        function robustness = adversarial_robustness_test_for_one_input(neural_net, input_vec, desired_output_val, attack_pos, perturbance_bound)
            %   @neural_net: the neural_net needs to test the robustness,
            %   currently only supports FFNN
            %   @input_vec: is the input vector corresponding to one image of one digit
            %   @desired_output_val: the range of the output that the neural network
            %   should produce
            %   @attack_pos: an array of attacked pixels positions
            %   @perturbance_bound: bound of perturbance caused by attacks

            %   @robustness:       = true -> the neural network produces correct output
            %                      = false -> the neural network produces incorrect
            %                      output (outside of the desired range)


            

        end

        
    end
    
end

