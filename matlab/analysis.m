format long
data = load("data/time_data.mat");
alice_times = data.alice_times;
bob_times = data.bob_times;

Ns = data.Ns;
Nb = data.Nb;

total_time_sec = data.label;
true_offset = data.offset;

disp(['Total protocol time (seconds): ', num2str(total_time_sec)]);
disp(['True offset: ', num2str(true_offset)]);


% Units
% 1 = micro sec
% 0.01 = nano sec
% 0.00001 = pico sec



N = 2^26;

total_time = double(total_time_sec * 1e6);


time_res = total_time/N;
disp(['Time resolution: ', num2str(time_res)])

analog_alice = to_analog(alice_times, N, total_time);
analog_bob = to_analog(bob_times, N, total_time);


lag_range = 5000;
lag_range_time = lag_range * time_res;
[cc, lags] = xcorr(analog_alice, analog_bob, lag_range);


plot(lags*time_res, cc);
xlabel('\tau (micro sec)');
ylabel('C(\tau)');
% ylim([ 1800, 2450]);

[peak, indx] = max(cc);
offset = lags(indx)*time_res;

disp(['Offset: ', num2str(offset)]);

mean_cc = mean(cc);
std_cc = std(cc);
S = (peak - mean_cc) / std_cc;
% S = (cc - mean_cc) / mean((cc - mean_cc).^2);

S_t = Ns /sqrt(double(Nb)/double(N));

p_err = N/2 * (1- erf(double(S)/sqrt(2)));

disp(['Coincidence frac: ', num2str(double(Ns)/double(Nb))]);
disp(['S: ', num2str(S)]);
% disp(['S_theory: ', num2str(S_t)]);
disp(['error probability: ', num2str(p_err)]);

c = 299792458;
disp(['Desired time accuracy: ', num2str((.1/c)*(1e6))]);
disp(['Accuracy: ', num2str(abs(true_offset+offset))])

% 
% bob_times = bob_times + offset;
% coincidences = [];
% for i=1:length(alice_times)
%     ti = alice_times(i);
%     bi = find(abs(bob_times - ti)<time_res);
%     if length(bi) == 1
%         coincidences = [coincidences, bob_times(bi(1))];
%     elseif length(bi) > 1
%         disp(['No single out possible for time', num2str(ti)]);
%     end
% end

