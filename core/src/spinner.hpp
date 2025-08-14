#ifndef LAHUTA_SPINNER_HPP
#define LAHUTA_SPINNER_HPP

#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>

// clang-format off

namespace lahuta {
namespace indicators {

struct SpinnerDesign {
    std::vector<std::string> frames;
};

inline std::unordered_map<std::string, SpinnerDesign> SpinnersMap = {
    {"none", {{""}}},
    {"dots", {{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"}}},
    {"ball", {{"( ●    )", "(  ●   )", "(   ●  )", "(    ● )", "(     ●)", "(    ● )", "(   ●  )", "(  ●   )", "( ●    )", "(●     )"}}},
    /*{"simple_dots", {{".  ", ".. ", "...", "   "}}},*/
    /*{"simple_dots_scrolling", {{".  ", ".. ", "...", " ..", "  .", "   "}}},*/
    /*{"star", {{"✶", "✸", "✹", "✺", "✹", "✷"}}},*/
    {"moon", {{"🌑 ", "🌒 ", "🌓 ", "🌔 ", "🌕 ", "🌖 ", "🌗 ", "🌘 "}}},
    {"hearts", {{"💛 ", "💙 ", "💜 ", "💚 ", "❤️ "}}},
    {"point", {{"∙∙∙", "●∙∙", "∙●∙", "∙∙●", "∙∙∙"}}},
    {"arrow2", {{"▹▹▹▹▹", "▸▹▹▹▹", "▹▸▹▹▹", "▹▹▸▹▹", "▹▹▹▸▹", "▹▹▹▹▸"}}},
};

class MinimalProgressSpinner {
public:
    MinimalProgressSpinner(const std::string& postfix_text, size_t max_progress)
        : prefix_text_(postfix_text), max_progress_(max_progress), tick_count_(0), currennt_state(0) {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, SpinnersMap.size() - 1);
        auto it = SpinnersMap.begin();
        std::advance(it, dist(gen));
        state = it->second.frames;
    }

    void tick() {
        std::lock_guard<std::mutex> lock(mutex_);
        tick_count_++;
        currennt_state = (currennt_state + 1) % state.size();
        print_progress();
    }

    void set_postfix_text(const std::string& postfix_text) { prefix_text_ = postfix_text; }

    void print_progress() {
        float percent_finished = 100.0f * tick_count_ / max_progress_;
        std::cout << "\r  " << state[currennt_state]
                  << " " << prefix_text_
                  << " " << std::fixed << std::setprecision(1) << percent_finished << "%"
                  << " (" << tick_count_ << "/" << max_progress_ << ")"
                  << std::flush;
    }

    void set_max_progress(size_t max_progress) { max_progress_ = max_progress; }
    std::mutex& get_mutex() { return mutex_; }
    void set_state(const std::string& state_name) { state = SpinnersMap[state_name].frames; }

private:
    std::vector<std::string> state;
    std::string prefix_text_;
    size_t max_progress_;
    size_t tick_count_;
    size_t currennt_state;
    std::mutex mutex_;
};


class SpinnerAwareSink : public spdlog::sinks::sink {
public:
    SpinnerAwareSink(std::shared_ptr<spdlog::sinks::sink> sink, MinimalProgressSpinner* spinner, std::mutex &spinner_mutex)
      : sink_(std::move(sink)), spinner_(spinner), spinner_mutex_(spinner_mutex) {}

    void log(const spdlog::details::log_msg &msg) override {
        if (spinner_) {
            std::lock_guard<std::mutex> lock(spinner_mutex_);   // lock the spinner: avoid log and spinner output mixing
            std::cout << "\r\033[K";                            // clears the current spinner line
            sink_->log(msg);                                    // pass the log message to the wrapped sink
            spinner_->print_progress();                         // reprint the spinner
        } else {
            sink_->log(msg);
        }
    }

    void flush() override { sink_->flush(); }
    void set_pattern(const std::string &pattern) override { sink_->set_pattern(pattern); }
    void set_formatter(std::unique_ptr<spdlog::formatter> sink_formatter) override {
        sink_->set_formatter(std::move(sink_formatter));
    }

private:
    std::shared_ptr<spdlog::sinks::sink> sink_;
    MinimalProgressSpinner* spinner_;
    std::mutex &spinner_mutex_;
};

} // namespace indicators
} // namespace lahuta

#endif // LAHUTA_SPINNER_HPP
