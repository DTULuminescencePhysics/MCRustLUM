use std::{error::Error, path::Path};

fn main() {
    monte_carlo_run();
   
}

fn monte_carlo_run() -> Result<(), Box<dyn Error>> {
    let inputs = if Path::new("input.toml").exists() {
        io::read_inputs("input.toml")?
    } else {
        io::default_inputs()
    };

    

    Ok(())
}
