#pragma once
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <csignal>

#ifdef __INTELLISENSE__
#define __FILE_NAME__ __FILE__
#endif

namespace neat {

class Logger {
public:
    static Logger& instance() {
        static Logger instance;
        return instance;
    }

    std::shared_ptr<spdlog::logger> get(const std::string& name) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = loggers_.find(name);
        if (it != loggers_.end()) {
            return it->second;
        }
        auto logger = create_logger(name);
        loggers_[name] = logger;
        return logger;
    }

    void set_level(const std::string& name, spdlog::level::level_enum level) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (auto it = loggers_.find(name); it != loggers_.end()) {
            it->second->set_level(level);
        }
    }

    void set_global_level(spdlog::level::level_enum level) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& [_, logger] : loggers_) {
            logger->set_level(level);
        }
    }

private:
    Logger() {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
            "logs/neat.log", 10 * 1024 * 1024, 5);

        std::string pattern = "%^%L%C%m%d|%H%M| %v%$";
		
		console_sink->set_pattern(pattern);
		file_sink->set_pattern(pattern);

        sinks_ = {console_sink, file_sink};
    }

    std::shared_ptr<spdlog::logger> create_logger(const std::string& name) {
        auto logger = std::make_shared<spdlog::logger>(name, sinks_.begin(), sinks_.end());
        logger->set_level(spdlog::level::trace);
        return logger;
    }

    std::vector<spdlog::sink_ptr> sinks_;
    std::unordered_map<std::string, std::shared_ptr<spdlog::logger>> loggers_;
    std::mutex mutex_;
};

} // namespace neat

#define LOGGER(name) neat::Logger::instance().get(name)
#define LOG_TRACE(...) do { if (logger) logger->trace("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_DEBUG(...) do { if (logger) logger->debug("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_INFO(...) do { if (logger) logger->info("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_WARN(...) do { if (logger) logger->warn("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_ERROR(...) do { if (logger) logger->error("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_CRITICAL(...) do { if (logger) logger->critical("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); } while(0)
#define LOG_FATAL(...) do { if (logger) { logger->critical("{:<{}} | {:>{}}", __VA_ARGS__, 186, fmt::format("{}:{}", __FILE_NAME__, __LINE__), 30); std::raise(SIGTERM); } } while(0)