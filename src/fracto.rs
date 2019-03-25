//! Generate a fracto image given the cycle by cycle crack growth
//! history.
//!
//! This currently creates a pseudo image where the color of each band
//! is a function of the Kmax. But in the future it should be a
//! rendered fracture surface based on the angle of growth.  The image
//! is output in svg format where each cycle is represented as a line
//! in the image. Because of this, it may be best to restrict the
//! number of cycles by limiting the starting and ending crack sizes.
//! The start of each block is indicated by drawing a red line through
//! the image. A periodic white scale bar is drawn at the bottom of
//! the image. The scale of the bar and the size of the image can be
//! changed.

use grow::History;
use std::str::FromStr;

// create the image using SVG
use svg::Document;
use svg::node::element;
use svg;

/// Data for generation of a reconstructed image.
#[derive(Debug, Clone)]
pub struct ImageData {
    /// Name of pseudo image file.
    pub file: String,
    /// Size of scale bar for image.
    pub barlength: f64,
    pub xsize: u32,
    pub ysize: u32,
    pub image: ImageType,
}

/// Image generation types.
#[derive(Debug, Clone)]
pub enum ImageType {
    Sem,
    Optical,
}

/// Implement the trait to get `ImageType` from a string.
impl FromStr for ImageType {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Sem" => Ok(ImageType::Sem),
            "Optical" => Ok(ImageType::Optical),
            _ => Err("no match"),
        }
    }
}

/// Write an image to a file.
// pub fn write(filename: &str, image: ImageBuffer<Luma<u8>, Vec<u8>>) {
//     let ref mut fout = File::create(&Path::new(filename)).unwrap();
 
//     // We must indicate the image's color type and what format to save as.
//     let _ = ImageLuma8(image).save(fout, PNG);
// }

// color models return a value  0.0 to 1.0
// simple color model based on the r ratio
//fn r_color_model(cycle: &History) -> f64 {
    // cap the r ratio
//    let r: f64 = if cycle.r < 0.0 {0.0} else {cycle.r};
//    r.powf(0.4)    
//}

/// Determines the color of a marker band.
/// Essentially the angle of the fracture plane to the average crack plane. Dark regions reflect the light away.
/// This model indicate the color is based on Kmax.
fn kmax_color_model(cycle: &History, cycle_max: &History) -> f64 {
    // cap the r ratio
    (cycle.k[0] / cycle_max.k[0]).max(0.0).powf(0.6).min(1.0)
}

/// Make a pseudo image of the fracture surface.
///
/// Fit the entire crack growth history into the image frame.
pub fn write_svg_pseudo_image(history: &[History], frame: &ImageData, filename: &str) {
    let mut document = Document::new().set("viewBox", (0, 0, frame.ysize, frame.xsize));

    let lightest = 1.0;
    println!("frame: {:?}", frame);
    // calculate the maximum cycle in the history
    //    let history_max = grow::cumulative_history(history);
    println!("***Tempory hack in image***");
    println!("cycles in history {}", history.len());
    let history_max = history;

    // check that the last point is the maximum
    let a_max = history[history.len() - 1].crack.a[0];
    let a_init = history[0].crack.a[0];
    println!("history {}, a_max {}", history.len(), a_max);

    let mut a_prev = a_init;
    println!("starting crack size is {}", a_prev);
    let mut delta_a = 0.0;
    let mut sum_color = 0.0;
    let mut ave_color;
    let mut y_prev = 0.0;
    let nbar = 10i32; // put this number of micron bars along the image

    //average over the length of a pixel
    let pixel_size = (a_max - a_init) / frame.ysize as f64; // size of a pixel

    let mut cblock = history[0].block.floor() + 1.0;
    let bar_pixels = (frame.barlength / pixel_size) as u32;

    // cover the canvas in black
    let path = element::Rectangle::new()
        .set("fill", "black")
        .set("stroke", 0)
        .set("y", 0)
        .set("x", 0)
        .set("width", frame.ysize)
        .set("height", frame.xsize);

    document = document.add(path);

    // plot each cycle in the history
    for (cycle, cycle_max) in history.iter().zip(history_max.iter()) {
        //        println!("{:?}", cycle);
        let striation = cycle.crack.a[0] - a_prev; // the width of the current striation
        a_prev = cycle.crack.a[0];

        // keep adding to the increment of growth until we have a single pixel big enough to display/
        if delta_a < pixel_size {
            delta_a += striation;
            sum_color += kmax_color_model(cycle, cycle_max) * striation; // weight the color according to the striation
            continue;
        }

        ave_color = sum_color / delta_a;
        sum_color = 0.0;
        delta_a = 0.0;

        let y_cur = frame.ysize as f64 * ((cycle.crack.a[0] - a_init) / (a_max - a_init));
        let newblock = if cycle.block > cblock {
            cblock = cycle.block.floor() + 1.0;
            true
        } else {
            false
        };

        // color in all the pixels from the previous pixel to the current pixel position
        let color = (lightest - ave_color).max(0.0);
        let shade = (color * 16.0) as u8;
        let code = format!("#{:x}{:x}{:x}", shade, shade, shade);
        let path = element::Rectangle::new()
            .set("fill", code)
            .set("stroke", 0)
            .set("y", 0)
            .set("x", y_prev)
            .set("width", (y_cur - y_prev))
            .set("height", frame.xsize - 50);

        document = document.add(path);

        // draw a line at the position of the start of the block
        if newblock {
            let start = element::Rectangle::new()
                .set("fill", "red")
                .set("y", 0.0)
                .set("x", y_prev)
                .set("height", frame.xsize)
                .set("width", 2.0);

            document = document.add(start);
        }

        y_prev = y_cur;
    }

    // add a micron bar every nbar spacings
    for b in 0..nbar {
        let bar_pos = b * frame.ysize as i32 / nbar;
        let bar = element::Rectangle::new()
            .set("fill", "white")
            .set("stroke", 0)
            .set("y", frame.xsize - 40)
            .set("x", bar_pos)
            .set("width", bar_pixels)
            .set("height", 10);

        document = document.add(bar);
    }

    // set the background to white
    // let path = element::Rectangle::new()
    //     .set("fill", "white")
    //     .set("stroke", 0)
    //     .set("y", frame.xsize - 2)
    //     .set("x", 0)
    //     .set("width", frame.ysize)
    //     .set("height", 2);

    // document = document.add(path);

    svg::save(filename, &document).unwrap();
}

// put a micron bar on the image
// fn pixel_micron_bar(imgx: u32, ypos: u32, imbuf: &mut ImageBuffer<Luma<u8>, Vec<u8>>, npixels: u32) {
//     // scale bar sizes
//     let bar_pwidth = 10u32;
//     let bar_color = Luma([255 as u8]);

//     // put micron along the length of the image
//     for z in 0u32..npixels {
//         for x in imgx - 5 - bar_pwidth..imgx - 5 {
//             imbuf.put_pixel(ypos + z, x, bar_color);
//         }
//     }
// }

// #[test]
// fn test_make_pseudo_image() {
//     let mut history = Vec::new();

//     for i in 0..10000 {
//         history.push((i as f64 / 100000.0,
//                       i as f64 / 100000.0,
//                       i as f64 / 100000.0,
//                       i as f64 / 10000.0));
//     }
//     let image = ImageData {file: "test.png".to_string(), barlength: 5e-6, xsize:300, ysize: 800};
//     make_pseudo_image(history, image, 5e-6);
// }
