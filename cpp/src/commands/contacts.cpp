#include "commands/contacts.hpp"
#include "cli/arg_validation.hpp"
#include "logging.hpp"

#include "db/db.hpp"
#include "entities/formatter.hpp"
#include "io/collector.hpp"
#include "io/file_backend.hpp"
#include "io/file_spill_policy.hpp"
#include "io/gzip_file_spill_policy.hpp"
#include "models/factory.hpp"
#include "pipeline/dsl.hpp"
#include "serialization/formats.hpp"
#include "tasks/contacts_task.hpp"
#include "tasks/model_db_writer.hpp"

#define STB_SPRINTF_IMPLEMENTATION
#include "gemmi/third_party/stb_sprintf.h"

// clang-format off
namespace lahuta::cli {

using namespace lahuta::pipeline;

using ContactsRes = tasks::ContactsTask::result_type;
using ptrT = std::shared_ptr<const ContactsRes>;
using Source = std::variant<sources::DirectorySource, sources::VectorSource, sources::FileListSource>;

static void initialize_factories(int num_threads) {
  lahuta::InfoPoolFactory::initialize(num_threads);
  lahuta::BondPoolFactory::initialize(num_threads);
  lahuta::AtomPoolFactory::initialize(num_threads);
}

struct ContactsOptions {
  enum class SourceMode { Directory, Vector, FileList, Database };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::string extension = ".cif";
  bool recursive = true;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;

  tasks::ContactProvider provider = tasks::ContactProvider::MolStar;
  InteractionType interaction_type = InteractionType::None;

  bool want_json   = false;
  bool want_text   = false;
  bool want_log    = false;
  bool no_compress = false;

  int threads = 8;
  size_t batch_size = 200;
};


Source pick_source(const ContactsOptions& cli) {
  switch (cli.source_mode) {
    case ContactsOptions::SourceMode::Directory: return sources::DirectorySource{cli.directory_path, cli.extension, cli.recursive, cli.batch_size};
    case ContactsOptions::SourceMode::Vector:    return sources::VectorSource   {cli.file_vector};
    case ContactsOptions::SourceMode::FileList:  return sources::FileListSource {cli.file_list_path};
    case ContactsOptions::SourceMode::Database:  break; // Handled separately
  }
  throw std::logic_error("Invalid source mode");
}

auto build_compute_stage(tasks::ContactProvider provider, InteractionType interaction_type) {
  tasks::ContactsTask task{provider, interaction_type};
  return stage(dsl::thread_safe, [task](std::string path, IEmitter<ptrT>& out) mutable {
    auto result = task(std::move(path));
    out.emit(std::make_shared<const ContactsRes>(std::move(result)));
  });
}

// Database deserialize stage
template<class Rec, class FormatTag>
auto make_database_deserialize_stage(lahuta::LMDBReader& rdr) {
    using Ptr = std::shared_ptr<const Rec>;
    return stage(dsl::thread_safe,
        [&rdr](std::string&& key, IEmitter<Ptr>& out)
        {
            std::string_view raw;
            if (!rdr.fetch(key, raw)) {
                 lahuta::Logger::get_logger()->error("missing key {}", key);
                return;
            }
            Rec obj = serialization::Serializer<FormatTag,Rec>
                          ::deserialize(raw.data(), raw.size());
            out.emit(std::make_shared<const Rec>(std::move(obj)));
        });
}

// Contacts computation stage that takes deserialized model data and uses ContactsTask
auto build_contacts_from_model_stage(tasks::ContactProvider provider, InteractionType interaction_type) {
  using WriterRes = tasks::ModelWriteTask::result_type;
  using InputPtr = std::shared_ptr<const WriterRes>;

  tasks::ContactsTask task{provider, interaction_type};
  return stage(dsl::thread_safe, [task](InputPtr input, IEmitter<ptrT>& out) {
    auto result = task(*input);
    out.emit(std::make_shared<const ContactsRes>(std::move(result)));
  });
}

// Output multiplexer components
class ICollectorSink {
public:
  virtual ~ICollectorSink() = default;
  virtual void emit(ptrT) = 0;
};

template<class CollectorT>
class CollectorSink final : public ICollectorSink {
public:
  template<class... Args>
  explicit CollectorSink(Args&&... args) : col_(std::forward<Args>(args)...) {}
  void emit(ptrT p) override { col_.emit(std::move(p)); }
  ~CollectorSink() override  { col_.finish(); }

private:
  CollectorT col_;
};

enum class OutFmt { Json, Text, Log };
struct OutputOption {
  OutFmt fmt;
  bool   compress;
};

template<OutFmt Fmt, bool Compress> struct CollectorTraits;

template<> struct CollectorTraits<OutFmt::Json, false> {
  using type = Collector<ptrT, FileSpillPolicy<fmt::json, ptrT, FileBackend, std::string>>;
  static constexpr std::string_view ext = ".json";
};
template<> struct CollectorTraits<OutFmt::Json, true>  {
  using type = Collector<ptrT, GzipFileSpillPolicy<fmt::json, ptrT>>;
  static constexpr std::string_view ext = ".json.gz";
};
template<> struct CollectorTraits<OutFmt::Text, false> {
  using type = Collector<ptrT, FileSpillPolicy<fmt::text, ptrT, FileBackend, std::string>>;
  static constexpr std::string_view ext = ".txt";
};
template<> struct CollectorTraits<OutFmt::Text, true>  {
  using type = Collector<ptrT, GzipFileSpillPolicy<fmt::text, ptrT>>;
  static constexpr std::string_view ext = ".txt.gz";
};

namespace {
  std::string format_string(const char* format, ...) {
    char buffer[2048];
    va_list args;
    va_start(args, format);
    stbsp_vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);
    return std::string(buffer);
  }

  std::string extract_filename(const std::string& filepath) {
    size_t pos = filepath.find_last_of("/\\");
    return (pos == std::string::npos) ? filepath : filepath.substr(pos + 1);
  }
}

class LoggingSink final : public ICollectorSink {
public:
  void emit(ptrT p) override {
    if (!p) return;

    const auto& result = *p;
    std::string filename = extract_filename(result.file_path);

    std::string output;
    std::string contact_type_str = (result.contact_type == InteractionType::None)
      ? "All"
      : interaction_type_to_string(result.contact_type);

    output += format_string("Contact Results for %s", filename.c_str()) + "\n";
    output += format_string("Provider: %s, Type: %s, Success: %s",
                           result.provider == tasks::ContactProvider::Arpeggio ? "Arpeggio" : "MolStar",
                           contact_type_str.c_str(),
                           result.success ? "Yes" : "No") + "\n";

    if (!result.success) {
      output += format_string("%s: Processing failed", filename.c_str()) + "\n";
    } else if (result.contacts.empty()) {
      output += format_string("%s: No contacts found", filename.c_str()) + "\n";
    } else {
      output += format_string("%s: Found %zu contacts:", filename.c_str(), result.contacts.size()) + "\n";
      for (const auto& contact : result.contacts) {
        if (result.topology) {
          std::string lhs_entity = ContactTableFormatter::format_entity_compact(*result.topology, contact.lhs);
          std::string rhs_entity = ContactTableFormatter::format_entity_compact(*result.topology, contact.rhs);
          output += format_string("%s:   %s <-> %s (distance: %.3f, type: %s)",
                                 filename.c_str(), lhs_entity.c_str(), rhs_entity.c_str(),
                                 contact.distance, interaction_type_to_string(contact.type).c_str()) + "\n";
        } else {
          output += format_string("%s:   %s <-> %s (distance: %.3f, type: %s)",
                                 filename.c_str(), contact.lhs.to_string().c_str(),
                                 contact.rhs.to_string().c_str(), contact.distance,
                                 interaction_type_to_string(contact.type).c_str()) + "\n";
        }
      }
    }
    output += "\n";

    std::cout << output << std::flush;
  }

  ~LoggingSink() override = default;
};

template<> struct CollectorTraits<OutFmt::Log, false> {
  using type = LoggingSink;
  static constexpr std::string_view ext = "";
};
template<> struct CollectorTraits<OutFmt::Log, true> {
  using type = LoggingSink;
  static constexpr std::string_view ext = "";
};

class OutputMultiplexer : public IEmitter<ptrT> {
  std::vector<std::unique_ptr<ICollectorSink>> sinks_;
public:
  explicit OutputMultiplexer(std::vector<OutputOption> opts, std::size_t batch, std::string_view stem) {
    if (opts.size() == 0) throw std::runtime_error("at least one output format required");

    for (std::size_t i = 0; i < opts.size(); ++i) {
      auto o = opts[i];
      switch (o.fmt) {
      case OutFmt::Json:
        (o.compress)
          ? add<OutFmt::Json, true >(batch, stem)
          : add<OutFmt::Json, false>(batch, stem);
        break;
      case OutFmt::Text:
        (o.compress)
          ? add<OutFmt::Text, true >(batch, stem)
          : add<OutFmt::Text, false>(batch, stem);
        break;
      case OutFmt::Log:
        // Compression doesn't apply to logging output
        add<OutFmt::Log, false>(batch, stem);
        break;
      }
    }
  }

  void emit(ptrT p) override {
    for (auto& s : sinks_) s->emit(p);
  }

private:
  template<OutFmt Fmt, bool Compress>
  void add(std::size_t batch, std::string_view stem) {
    using ColT = typename CollectorTraits<Fmt, Compress>::type;

    if constexpr (Fmt == OutFmt::Log) {
      // LoggingSink doesn't need file parameters
      sinks_.push_back(std::make_unique<LoggingSink>());
    } else {
      auto filename = std::string(stem).append(CollectorTraits<Fmt, Compress>::ext);
      sinks_.push_back(std::make_unique<CollectorSink<ColT>>(batch, filename));
    }
  }
};

template<typename StageT, typename SinkT>
void run_dispatch(Source& src_variant, StageT& stage, SinkT& sink, int threads) {
  std::visit([&](auto&& src) {
    auto pipeline = dsl::source_t{std::move(src)} | stage | dsl::collect(sink);
    run(pipeline, threads);
  }, src_variant);
}

namespace contacts_opts {
const option::Descriptor usage[] = {
  {ContactsOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta contacts [options]\n\n"
   "Compute inter-atomic contacts.\n\n"
   "Input Options (choose one):"},
  {ContactsOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {ContactsOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {ContactsOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated)."},
  {ContactsOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {ContactsOptionIndex::SourceDatabase, 0, "", "database", validate::Required,
   "  --database <path>            \tProcess structures from database."},
  {ContactsOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension for directory mode (default: .cif)."},
  {ContactsOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {0, 0, "", "", option::Arg::None,
   "\nCompute Options:"},
  {ContactsOptionIndex::Provider, 0, "p", "provider", validate::Provider,
   "  --provider, -p <provider>    \tContact provider: 'molstar' or 'arpeggio' (default: molstar)."},
  {ContactsOptionIndex::InteractionType, 0, "i", "interaction", validate::ContactType,
   "  --interaction, -i <type>     \tInteraction type: 'hbond', 'hydrophobic', 'ionic', etc.\n"},
  {0, 0, "", "", option::Arg::None,
   "\nOutput Options:"},
  {ContactsOptionIndex::OutputJson, 0, "", "json", option::Arg::None,
   "  --json                       \tOutput results in JSON format."},
  {ContactsOptionIndex::OutputText, 0, "", "text", option::Arg::None,
   "  --text                       \tOutput results in text format."},
  {ContactsOptionIndex::OutputLog, 0, "", "log", option::Arg::None,
   "  --log                        \tOutput results to standard output (logging)."},
  {ContactsOptionIndex::NoCompress, 0, "", "no-compress", option::Arg::None,
   "  --no-compress                \tDisable gzip compression (default: enabled)."},
  {0, 0, "", "", option::Arg::None,
   "\nRuntime Options:"},
  {ContactsOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {ContactsOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for processing (default: 200)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace contacts_opts

[[nodiscard]] std::unique_ptr<CliCommand> ContactsCommand::create() {
  return std::unique_ptr<CliCommand>(new ContactsCommand());
}

int ContactsCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, contacts_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, contacts_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[contacts_opts::ContactsOptionIndex::Help]) {
    option::printUsage(std::cout, contacts_opts::usage);
    return 0;
  }

  try {
    // parse CLI Arguments into ContactsOptions
    ContactsOptions cli;

    // parse source options
    int source_count = 0;
    if (options[contacts_opts::ContactsOptionIndex::SourceDirectory]) {
      cli.source_mode = ContactsOptions::SourceMode::Directory;
      cli.directory_path = options[contacts_opts::ContactsOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceVector]) {
      cli.source_mode = ContactsOptions::SourceMode::Vector;
      std::string files_str = options[contacts_opts::ContactsOptionIndex::SourceVector].arg;

      // parse comma-separated files
      std::stringstream ss(files_str);
      std::string file;
      while (std::getline(ss, file, ',')) {
        if (!file.empty()) {
          cli.file_vector.push_back(file);
        }
      }
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceFileList]) {
      cli.source_mode = ContactsOptions::SourceMode::FileList;
      cli.file_list_path = options[contacts_opts::ContactsOptionIndex::SourceFileList].arg;
      source_count++;
    }
    if (options[contacts_opts::ContactsOptionIndex::SourceDatabase]) {
      cli.source_mode = ContactsOptions::SourceMode::Database;
      cli.database_path = options[contacts_opts::ContactsOptionIndex::SourceDatabase].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --directory, --files, --file-list, or --database");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    // Parse other options
    if (options[contacts_opts::ContactsOptionIndex::Extension]) {
      cli.extension = options[contacts_opts::ContactsOptionIndex::Extension].arg;
    }

    cli.recursive = options[contacts_opts::ContactsOptionIndex::Recursive] ? true : false;

    // Parse provider
    if (options[contacts_opts::ContactsOptionIndex::Provider]) {
      std::string provider = options[contacts_opts::ContactsOptionIndex::Provider].arg;
      if (provider == "molstar") {
        cli.provider = tasks::ContactProvider::MolStar;
      } else if (provider == "arpeggio") {
        cli.provider = tasks::ContactProvider::Arpeggio;
      } else {
        Logger::get_logger()->error("Invalid provider '{}'. Must be 'molstar' or 'arpeggio'", provider);
        return 1;
      }
    }

    // Parse interaction type
    if (options[contacts_opts::ContactsOptionIndex::InteractionType]) {
      std::string interaction = options[contacts_opts::ContactsOptionIndex::InteractionType].arg;
      cli.interaction_type = get_interaction_type(interaction);
    }

    // Parse output options
    cli.want_json   = options[contacts_opts::ContactsOptionIndex::OutputJson] ? true : false;
    cli.want_text   = options[contacts_opts::ContactsOptionIndex::OutputText] ? true : false;
    cli.want_log    = options[contacts_opts::ContactsOptionIndex::OutputLog]  ? true : false;
    cli.no_compress = options[contacts_opts::ContactsOptionIndex::NoCompress] ? true : false;

    if (!cli.want_json && !cli.want_text && !cli.want_log) cli.want_json = true; // json if nothing specified

    // runtime options
    if (options[contacts_opts::ContactsOptionIndex::Threads]) {
      cli.threads = std::stoi(options[contacts_opts::ContactsOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[contacts_opts::ContactsOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[contacts_opts::ContactsOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    // build and run pipeline
    auto compute_stage = build_compute_stage(cli.provider, cli.interaction_type);

    // output multiplexer
    std::vector<OutputOption> option_list;
    if (cli.want_json) { option_list.push_back(OutputOption{OutFmt::Json, !cli.no_compress}); }
    if (cli.want_text) { option_list.push_back(OutputOption{OutFmt::Text, !cli.no_compress}); }
    if (cli.want_log)  { option_list.push_back(OutputOption{OutFmt::Log, false}); } // Compression n/a for logging
    OutputMultiplexer output_emitter{std::move(option_list), cli.batch_size, "contacts"};

    if (cli.source_mode == ContactsOptions::SourceMode::Database) {
      initialize_factories(cli.threads);

      // special handling for database sources
      LMDBDatabase database(cli.database_path);
      LMDBReader reader(database.get_env(), database.get_dbi());

      // we create the pipeline
      auto key_source         = dsl::db_keys(database, cli.batch_size);
      auto deserialize_stage  = make_database_deserialize_stage<tasks::ModelWriteTask::result_type, fmt::binary>(reader);
      auto compute_stage      = build_contacts_from_model_stage(cli.provider, cli.interaction_type);

      auto pipeline = key_source | deserialize_stage | compute_stage | output_emitter;
      dsl::run(pipeline, cli.threads);
    } else {
      // file-based sources
      Source source_variant = pick_source(cli);
      run_dispatch(source_variant, compute_stage, output_emitter, cli.threads);
    }

    Logger::get_logger()->info("Contact computation completed successfully!");
    return 0;

  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
