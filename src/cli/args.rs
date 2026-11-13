//! Options that several commands share. Each set is flattened into every command that takes
//! it, so that an option has the same name, short form and help text everywhere; an option
//! that only one solver has is declared by that solver's command.

use super::Error;
use omelettehedron::cube::MoveSet;
use std::path::PathBuf;
use std::time::Duration;

/// The positions to work on (`solve`, `show`).
#[derive(clap::Args, Clone, Debug)]
#[command(next_help_heading = "Input")]
pub struct InputArgs {
    /// A scramble in FTO notation (U R F L B BR BL D, ' for anticlockwise; a w after the face
    /// for a wide turn, an s for a slice turn), or a pattern: the 72 sticker colours, as
    /// `omelettehedron show --help` describes. Without one (and without --random or
    /// --patterns): one per line on stdin, `#` starting a comment
    #[arg(value_name = "SCRAMBLE|PATTERN")]
    pub position: Vec<String>,
    /// Each scramble given is a solution: use the position it solves (the scramble inverted,
    /// applied to the solved puzzle) instead of the position it makes
    #[arg(long, conflicts_with_all = ["random", "patterns"])]
    pub inverse: bool,
    /// Read patterns from stdin, each on one line or several (lines are joined until they hold
    /// 72 stickers; a blank line ends one early); `#` starts a comment. A pattern with a
    /// mistake is reported and the others are still used
    #[arg(long, conflicts_with_all = ["position", "random"])]
    pub patterns: bool,
    /// Use this many seeded random positions
    #[arg(long, value_name = "N", conflicts_with = "position")]
    pub random: Option<usize>,
    /// Seed of the random positions (--random)
    #[arg(long, default_value_t = 1, requires = "random")]
    pub seed: u64,
}

/// `--usual`, for the scramble printed after the solution of a pattern (`solve`).
#[derive(clap::Args, Clone, Debug)]
#[command(next_help_heading = "Input")]
pub struct UsualArgs {
    /// For positions given as patterns: your solved puzzle as you usually hold it to scramble,
    /// as the colours of U R F L B BL D BR (8 characters, the letters of the pattern). The
    /// scramble printed after the solution then also starts from that hold
    #[arg(long, value_name = "COLOURS")]
    pub usual: Option<String>,
}

/// When a solve stops (`solve`, `bench`); without either, it runs until its solution is
/// proven optimal.
#[derive(clap::Args, Clone, Debug)]
#[command(next_help_heading = "Limits")]
pub struct LimitArgs {
    /// Stop after this many seconds: the best solution found by then, and the lower bound
    /// proven
    #[arg(short = 't', long, value_name = "SECONDS")]
    pub time: Option<f64>,
    /// Stop at the first solution of at most this many moves
    #[arg(short = 'l', long, value_name = "MOVES")]
    pub max_length: Option<usize>,
}

impl LimitArgs {
    /// The time limit, if any.
    pub fn duration(&self) -> Result<Option<Duration>, Error> {
        seconds(self.time, "--time")
    }
}

/// `secs` as a duration (`None` stays `None`); errors name the option `name`.
pub fn seconds(secs: Option<f64>, name: &str) -> Result<Option<Duration>, Error> {
    secs.map(Duration::try_from_secs_f64)
        .transpose()
        .map_err(|e| format!("{name}: {e}").into())
}

/// Where the tables are, which moves solutions use, and how tables may be generated (every
/// command that loads tables). The tier is each command's own (each solver has its tiers).
#[derive(clap::Args, Clone, Debug)]
#[command(next_help_heading = "Tables")]
pub struct TableArgs {
    /// Directory of the pruning tables (default: $FTO_TABLES or ./tables). Missing tables of
    /// the tier are generated there
    #[allow(clippy::doc_markdown)] // shown verbatim by --help
    #[arg(long, value_name = "DIR")]
    pub tables: Option<PathBuf>,
    /// Move families allowed in solutions besides face turns, comma-separated; every move
    /// counts as one. slice: middle-layer turns (own tables, `*-slice.prun`). wide: two-layer
    /// turns; they never shorten a solution (a wide turn is the opposite face turn plus a
    /// rotation), so only the output changes: wide turns of the back faces instead of
    /// front-face turns, with the front corner held still
    #[arg(long, value_enum, value_delimiter = ',')]
    pub moves: Vec<Family>,
    /// Generate missing tables even if the disk or memory check says they do not fit
    #[arg(long)]
    pub force: bool,
    /// Memory (GB) to leave to the rest of the machine: tables are generated only if this much
    /// stays free beyond them, and the generation of the optimal solver's large tier stops
    /// (writing nothing) if free memory falls below it
    #[arg(long, value_name = "GB", default_value_t = 0.0)]
    pub min_free_memory: f64,
}

/// A move family for `--moves`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Family {
    /// two-layer turns, e.g. `Uw`
    Wide,
    /// middle-layer turns, e.g. `Us` (= `Ds'`)
    Slice,
}

impl TableArgs {
    /// The move set of `--moves`.
    pub fn move_set(&self) -> MoveSet {
        MoveSet {
            wide: self.moves.contains(&Family::Wide),
            slice: self.moves.contains(&Family::Slice),
        }
    }

    /// Set the checks before generating: `--force` and `--min-free-memory`.
    pub fn apply_checks(&self) {
        omelettehedron::preflight::set_force(self.force);
        // (an `as` cast saturates: negative or NaN -> 0, too large -> u64::MAX)
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let floor = (self.min_free_memory * 1e9) as u64;
        omelettehedron::preflight::set_memory_floor(floor);
    }
}

/// What a solve prints besides its result (`solve`).
#[derive(clap::Args, Clone, Debug)]
#[command(next_help_heading = "Output")]
pub struct OutputArgs {
    /// Print only each position's final solution (one line each), and nothing on stderr but
    /// errors
    #[arg(short, long, conflicts_with_all = ["verbose", "progress"])]
    pub quiet: bool,
    /// Also print each position's name, the solver's statistics (the optimal solver's: also
    /// its bounds at the start and each iteration) and, for a scramble, whether the solution
    /// also solves the pure puzzle (every centre its own piece)
    #[arg(short, long)]
    pub verbose: bool,
    /// Print what the search is doing to stderr while it runs: each iteration as it ends
    /// (depth, nodes, nodes/s, best length, lower bound proven) and, during a long iteration, a
    /// line every SECONDS (default 1; 0: none) with the nodes so far and the tasks started
    #[arg(
        long,
        value_name = "SECONDS",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "1"
    )]
    pub progress: Option<f64>,
}

impl OutputArgs {
    /// The heartbeat interval of `--progress`, if progress is on.
    pub fn every(&self) -> Result<Option<Duration>, Error> {
        seconds(self.progress, "--progress")
    }
}
