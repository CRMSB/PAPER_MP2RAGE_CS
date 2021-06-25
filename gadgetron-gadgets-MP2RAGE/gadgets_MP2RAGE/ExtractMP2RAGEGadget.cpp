/*
 * ExtractMagnitudeGadget.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */


#include <bitset>
#include <gadgetron/hoNDArray_math.h>
#include <unordered_map>

#include "ExtractMP2RAGEGadget.h"
#include <boost/math/constants/constants.hpp>

#include "hoNDArray_fileio.h"
#include <gadgetron/mri_core_def.h>


namespace Gadgetron {

    namespace {
        using IMTYPE = ISMRMRD::ISMRMRD_ImageTypes;

        const std::map<IMTYPE, size_t> series_offset{ { IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE, 0 },
            { IMTYPE::ISMRMRD_IMTYPE_REAL, 1000 }, { IMTYPE::ISMRMRD_IMTYPE_IMAG, 2000 },
            { IMTYPE::ISMRMRD_IMTYPE_PHASE, 3000 } };

        const std::map<int, IMTYPE> mask_to_imtype{ { 0, IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE },
            { 1, IMTYPE::ISMRMRD_IMTYPE_REAL }, { 2, IMTYPE::ISMRMRD_IMTYPE_IMAG },
            { 3, IMTYPE::ISMRMRD_IMTYPE_PHASE } };

        template <class FUNCTION>
        hoNDArray<float> extract(const hoNDArray<float>& data, FUNCTION&& extractor) {
            hoNDArray<float> output(data.dimensions());
            std::transform(data.begin(), data.end(), output.begin(), extractor);
            return output;
        }

        hoNDArray<float> extract_image(const hoNDArray<float>& data, IMTYPE imtype, float offset) {
            switch (imtype) {
            case ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE: return extract(data, [](auto& val) { return std::abs(val); });
            case ISMRMRD::ISMRMRD_IMTYPE_PHASE:
                return extract(data, [](auto& val) { return std::arg(val) + boost::math::constants::pi<float>(); });
            case ISMRMRD::ISMRMRD_IMTYPE_REAL: return extract(data, [&](auto& val) { return std::real(val) + offset; });
            case ISMRMRD::ISMRMRD_IMTYPE_IMAG: return extract(data, [&](auto& val) { return std::imag(val) + offset; });
            default: throw std::runtime_error("Illegal image tpye encounted in extract_image");
            }
        }

    }

    void ExtractMP2RAGEGadget::process(Core::InputChannel<Core::Image<float>>& in, Core::OutputChannel& out) {
        GDEBUG("----------EXTRACT_MP2RAGE_PROCESS----------------\n");

        for (auto image : in) {
            int count =0;
            count = count +1;
            GDEBUG("----------Image in loop = %d----------------\n",count);
            const auto& head = std::get<ISMRMRD::ImageHeader>(image);
            const auto& data = std::get<hoNDArray<float>>(image);
            const auto& meta = std::get<2>(image);

            ISMRMRD::MetaContainer meta_copy = meta.value_or(ISMRMRD::MetaContainer());

            for (auto imtype : image_types) {
                auto data_copy       = extract_image(data, imtype, real_imag_offset);
                auto head_copy       = head;
                
                if (head.image_series_index == 5000)
                {
                  GDEBUG("------TI1 IMAGE------\n");
                  meta_copy.set(GADGETRON_IMAGECOMMENT,"Magn_TI1_CS");
                  meta_copy.set(GADGETRON_SEQUENCEDESCRIPTION, "Magn_TI1_CS");
                  meta_copy.set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
                }
                else if (head.image_series_index == 6000)
                {
                  GDEBUG("------TI2 IMAGE------\n");
                  meta_copy.set(GADGETRON_IMAGECOMMENT, "Magn_TI2_CS");
                  meta_copy.set(GADGETRON_SEQUENCEDESCRIPTION, "Magn_TI2_CS");
                  meta_copy.set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
                }
                else if (head.image_series_index ==7000)
                {
                  GDEBUG("------MP2RAGE IMAGE------\n");
                  meta_copy.set(GADGETRON_IMAGECOMMENT, "MP2RAGE");
                  meta_copy.set(GADGETRON_SEQUENCEDESCRIPTION, "MP2RAGE");
                  meta_copy.set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
                }
                else if(head.image_series_index ==8000)
                {
                  GDEBUG("------T1 map------\n");
                  meta_copy.set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);//
                  meta_copy.set(GADGETRON_IMAGECOMMENT, "T1_CS");
                  meta_copy.set(GADGETRON_SEQUENCEDESCRIPTION, "T1_CS");
                  //metaset(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
                  meta_copy.set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
                  //metaset(GADGETRON_IMAGE_COLORMAP,  "HotMetal.pal");
                  meta_copy.set(GADGETRON_IMAGE_COLORMAP,  "Perfusion.pal");//
                }
                else
                {
                  GDEBUG("------Undefined------\n");
                  meta_copy.set(GADGETRON_IMAGECOMMENT, "Undefined");
                  meta_copy.set(GADGETRON_SEQUENCEDESCRIPTION, "Undefined");
                  meta_copy.set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
                }
                
                head_copy.data_type  = ISMRMRD::ISMRMRD_FLOAT;
                head_copy.image_type = imtype;
                head_copy.image_series_index += series_offset.at(imtype);
                out.push(head_copy, std::move(data_copy), meta_copy);
            }
        }
    }

    ExtractMP2RAGEGadget::ExtractMP2RAGEGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<Core::Image<float>>(context, props) {

        for (int i = 0; i < extract_mask.size(); i++) {
            if (extract_mask[i])
                image_types.insert(mask_to_imtype.at(i));
        }
        if (extract_magnitude)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE);
        if (extract_real)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_REAL);
        if (extract_imag)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_IMAG);
        if (extract_phase)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_PHASE);

        if (image_types.empty())
            throw std::runtime_error("ExctractGadget: No valid extract functions specified");
    }

    GADGETRON_GADGET_EXPORT(ExtractMP2RAGEGadget)

}
